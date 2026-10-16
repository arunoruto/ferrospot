# Code review: speed and correctness

Measured on 2026-10-16 (one machine, one run each, cProfile overhead
included; `pytest -m "not data"`: 31 passed). Nothing here is applied
yet; it is the list to work through, ranked by gain over risk.

## Where the time goes

A reduced calibration (one thickness, one seed, three states, 150 px
images, 13 simulated images plus 8 pilot runs) takes 3.8 s:

| step | share | note |
|---|---|---|
| `simulate_image`, thickness field | 47 % | one `gaussian_filter` with $\sigma = n/8$ on the 3× oversampled grid: a 360-tap kernel per axis |
| `detect` (13 calls) | 24 % | LoG 0.3 s, `half_max_radii` 0.42 s |
| `polydisc.sample` (11 calls) | 15 % | 0.05 s per saturated packing, 0.002 s at large gaps |
| `iron_path` (21 calls) | 8 % | 45 to 60 µs per sphere |

Single calls: `simulate_image` 0.29 s at 300 px (`iron_path` 0.26 s for
4400 spheres); `detect` DoG 0.20 s for 665 particles, of which
`half_max_radii` is 86 %; `detect` LoG with `split_touching` 0.94 s, of
which splitting is 56 % and includes a second full `half_max_radii`
pass only to get contrasts; `blob_log` on 1200 × 800 px with 20 scales
0.4 s on a real-looking image but 4.6 s on noise (pruning of local
maxima); `spacing_stats` + `neighbours` 8 ms (not a hot spot).

## Ranked changes

1. **Thickness field in `simulate_image`** (`simulate.py`, the
   `gaussian_filter(rng.normal(...), n/8, mode="wrap")`). 40 ms per call
   at 360 px, nearly half of a calibration. Options: (a) cache the field
   per `(shape, seed)`, since `match_noise` and `run_one` draw the same
   one again and again; (b) filter through `scipy.ndimage.fourier_gaussian`
   (same numbers up to rounding); (c) draw the field on the coarse grid
   and upsample, 20× less work but a different realisation, so the pinned
   regression numbers move. Do (a) or (b) first.
2. **`half_max_radii` loop** (`sizes.py`). 0.2 ms per particle:
   `np.indices`, `hypot`, `percentile`, `label`, `unravel_index` per
   window. Precompute one radial template per window size, `np.partition`
   instead of `percentile`, `ogrid` instead of `indices`: 30 to 40 %
   off, numbers unchanged. Preferred: a numba kernel over all particles
   (numba is already installed through polydisc; make it a core
   dependency or an optional accelerator with a NumPy fallback), which
   removes the Python loop entirely.
3. **`_split_all`** (`detect.py`). Drop the contrast pre-pass by taking
   the contrast from the final sizing, and call `h_maxima` only for blobs
   large enough to hold two particles (area above $1.5\pi r_\text{min}^2$).
   Up to 30 % of the LoG call. Medium risk: the pre-pass measures
   contrast on the normalised image, the final pass on the flattened one,
   so the split level moves slightly.
4. **`iron_path`** (`simulate.py`). `np.minimum`/`np.maximum` instead
   of `np.clip`, slice adds instead of `np.ix_` when a sphere does not
   wrap: 2 to 3× for free, bit-identical. Preferred: a numba kernel,
   about 50×, which also makes dense packings (2 to 4 s at 40 000 to
   80 000 spheres today) and finer oversampling cheap.
5. **`blob_log` scales**. 10 scales instead of 20 halves the time but
   changes which radii and blobs survive pruning; check against the
   regression panels before changing `NUM_SCALES`.
6. Not worth touching: packing cost in `_states` (parallelised by
   `--workers` already), the distance transform in `neighbours` (2.5 ms,
   three calls), the `overlay` circle loop.

## Correctness findings

Confirmed by experiment: `SIGMA_PER_RADIUS = 0.60` (the scale-normalised
LoG of a projected sphere peaks at $\sigma/R$ = 0.59 to 0.60 for $R$ = 3
to 10 px, with or without a 2 px blur); the Clark-Evans expectation
$0.5/\sqrt{N_A}$ (0.5000 measured) and the 3D Poisson constant 0.554
(0.5538 measured); `calibrate._match` cannot count one detection for two
spheres.

Half-maximum sizing against $R\sqrt3/2$, no noise:

| $R$ (px) | blur 0 | blur 2 px | blur 4 px |
|---|---|---|---|
| 2 | 0.98× | 0.98× | 1.49× |
| 3 | 1.00× | 1.00× | 1.09× |
| 5 | 0.87× | 0.87× | 0.87× |
| 10 | 0.91× | 0.91× | 0.91× |
| 20 | 0.93× | 0.92× | 0.92× |

Radii of resolved spheres read 7 to 13 % below the $\sqrt3/2$ outline
(pixel quantisation and the pre-smoothing), small blurred ones read up to
50 % large, as the module docstring says.

To fix, in order:

1. **Border correction is biased** (`stats.neighbours`, `counted = edge
   > nn`). Keeping a particle only when its own NN distance is shorter
   than its edge distance selects short distances: on a Poisson pattern
   the corrected mean NN reads 3.0 % low and the Clark-Evans ratio 0.972
   instead of 1 (uncorrected: 2.2 % high). Within the 10 % tolerance the
   calibration assumes, but systematic. Replace by a fixed guard distance
   (count particles farther than $d_\text{guard}$ from the edge, with
   neighbours searched among all) or the Hanisch estimator; the
   calibration must then use the same estimator on simulated images,
   which it does automatically through `spacing_stats`.
2. **Half-pixel offset.** `iron_path` samples pixel centres at
   $(i + 0.5)\,s$, `detect` reports $x = i\,s$: simulated positions and
   detections differ by a constant $-0.5$ px (median offset $-0.43,
   -0.37$ px measured, not a scale error). Spacings are unaffected;
   `overlay`, `_match` and any future labels are. Report
   $x = (i + 0.5)\,s$ in `detect` and mirror it in `_pixels`.
3. **`_pixels` truncation.** `int(x_nm / s)` lands one pixel low in about
   4 % of cases because `(i s)/s$ rounds below $i$; mask and depth lookups
   then use the neighbouring pixel. Use `floor(x / s + 1e-9)` or round
   after the half-pixel fix.
4. **Two centres 1 px apart** (LoG can produce them): the Voronoi cut
   halves the first particle and drops the second as NaN, silently. Merge
   detections closer than 1 px before sizing.
5. **Border particles** have their half-maximum region truncated by the
   image edge (radius 12 to 32 % low for a sphere cut by the edge) and
   are not flagged. Flag particles within $r_\text{max}$ of the edge or
   exclude them from the size statistics as the border correction does
   for spacings.
6. **Simulated pixel size** is off by up to 0.15 % when the oversampled
   width is not a multiple of `oversample` (12 % of the size/pixel
   combinations the calibration produces), with two fine pixels dropped at
   the seam. Choose the fine grid as a multiple of `oversample`.
7. **LoG and DoG thresholds are not interchangeable**: at the same
   `threshold` the DoG is 2.3× stricter (response 0.25 against 0.57 to
   0.63 for the same particle), and a DoG tuned to $r_\text{min}$
   responds to a particle of $3\,r_\text{min}$ with a fifth of that. The
   docs now say so; a scale-normalised DoG ($\times\sigma^2$ and the
   1/(k−1) factor) would make one threshold mean the same in both.
