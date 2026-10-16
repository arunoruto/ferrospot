# How the pipeline works

One image goes through five steps: load and scale, mask, detect, size,
measure. Every step is a plain function in its own module, so each can be
run, tested and replaced on its own. The equations below are the ones the
code implements; the symbols are the ones used in the source.

```{figure} _static/pipeline.svg
:alt: load → mask → detect → size → stats

Modules and what flows between them. Grey arrows are optional inputs.
```

## 1. Load and scale (`ferrospot.image`, `ferrospot.pdf`)

`load_image` returns a 2D `float64` array. Unsigned 8- and 16-bit files are
divided by their type maximum, so intensities are in $[0, 1]$ whatever the
file depth; float and 32-bit integer files keep their values. Colour is
converted to luminance, and only the first frame of a multi-frame TIFF is
read.

The pixel size $s$ (nm per pixel) comes from one of

* the microscope metadata (best, type it in);
* a scale bar measured by hand (`bar_nm / bar_px`);
* `measure_scale_bar`: threshold the search box at the midpoint of its
  intensity range, label connected components, take the widest one as the
  bar, so $s = L_\text{bar} / w_\text{px}$. A bar tilted by a degree or
  two still measures right because only the horizontal extent counts;
* for figures in PDFs, a vector bar measured in a page render at `dpi`
  and converted with the image's placement resolution (`pdf.image_ppi`):
  $s = L_\text{bar} / w_\text{render}\cdot \text{dpi}/\text{ppi}$.

## 2. Mask (`ferrospot.masks`)

`build_mask` returns the region to count in and, with a vacuum level,
the depth of every pixel. The mask is the intersection of

* an ROI box and the complement of excluded boxes (labels, insets, bars);
* an intensity window $[I_\text{lo}, I_\text{hi}]$ on the image smoothed
  with a Gaussian of $\sigma = $ `smooth_nm` $/ s$, with holes filled (the
  particles inside the window region) and eroded by `erode_nm`: the glass
  between dark vacuum and a bright mineral, for example;
* not the vacuum: the largest connected region below `vacuum_below`
  (dark field, HAADF) or above `vacuum_above` (bright field) in the
  smoothed image. Its Euclidean distance transform times $s$ is the
  depth map $z(x, y)$, and `depth_band_nm` keeps $z_0 \le z < z_1$.

Intensities here are in the loaded image's units, before any bright-field
inversion.

## 3. Detect (`ferrospot.detect`)

### Normalisation

`normalise` maps the 0.5 and 99.9 percentiles inside the mask to 0 and 1
(clipped), inverts bright-field images so particles are always bright,
and fills every pixel outside the mask with the value of the nearest pixel
inside (nearest-neighbour extrapolation through the distance transform).
The fill matters: without it a vacuum edge or a cut through a thickness
ramp is the strongest "blob" in the image. After this step `threshold`
means roughly the same thing across images, but it still sets the count:
try a few values and look at the overlay.

### LoG blobs (mixed sizes)

`detector="log"` runs `skimage.feature.blob_log` on the normalised image:
the scale-normalised Laplacian of Gaussian

$$
L(x, y; \sigma) = \sigma^2 \, \nabla^2 (G_\sigma * I)(x, y)
$$

over `NUM_SCALES` $= 20$ scales between $\sigma_\text{min}$ and
$\sigma_\text{max}$, and keeps local maxima of $-L$ in $(x, y, \sigma)$
above `threshold`. For a flat disc of radius $R$ the response peaks at
$\sigma = R / \sqrt 2$ {cite:p}`Lindeberg1998`. A projected sphere is not
a flat disc: its HAADF intensity goes with the chord length
$2\sqrt{R^2 - \rho^2}$, and for that profile the peak lies at

$$
\sigma^\star = 0.60\, R \qquad (\texttt{SIGMA\_PER\_RADIUS}),
$$

so `radius_nm=(r_min, r_max)` becomes $\sigma \in [0.6\, r_\text{min}/s,\;
0.6\, r_\text{max}/s]$ and a blob at scale $\sigma$ gets the rough radius
$\sigma / 0.6$. That radius is only a guess for the sizing step, not the
reported size.

`split_touching=True` looks inside every LoG blob for touching particles:
threshold the smoothed image at half the blob's contrast, take the
connected component under the centre, and if the distance transform of
that component has more than one $h$-maximum ($h = 1$ px) watershed it
into parts, each sized from the distance maximum as $R = d_\max /
(\sqrt 3 / 2)$ (see below for the $\sqrt 3 / 2$).

### DoG maxima (dense, small particles)

`detector="dog"` is one band-pass tuned to the smallest radius,

$$
B = G_{\sigma} * I - G_{1.6\,\sigma} * I, \qquad \sigma = 0.60\, r_\text{min}/s,
$$

and keeps local maxima of $B$ above `threshold` that are at least
$r_\text{min}/s$ pixels apart (`peak_local_max`). The ratio 1.6
(`DOG_RATIO`) makes the DoG the usual approximation of the LoG at that
scale, but the DoG is *not* scale-normalised, so a `threshold` of 0.02 in
DoG is not the same contrast as 0.02 in LoG. DoG separates dense particles
that the multi-scale LoG merges into one large blob; it cannot size them
(every detection starts with the guess $r_\text{min}$) and relies on the
sizing step for that.

Detections whose centre pixel is outside the mask are dropped.

## 4. Size (`ferrospot.sizes`)

Every particle, whichever detector found it, is sized the same way:
`half_max_radii` grows the connected region around the centre that is
brighter than halfway between the particle's peak and its local
background, and reports the equivalent-circle radius

$$
r = \sqrt{A_{1/2} / \pi}.
$$

Details that matter:

* the image is first flattened, $I - G_{2 r_\text{max}} * I$, removing
  thickness gradients and bright-field shading on scales well above the
  particles, then smoothed with $\sigma = \max(1, 0.25\, r_\text{min}/s)$
  so pixel noise does not fray the outline;
* the peak is the maximum within $r_\text{guess}/2$ of the centre; the
  background is the 10th percentile (`BG_PERCENTILE`) of the ring $1.5\,
  r_\text{guess} \le \rho \le 2.5\, r_\text{guess}$ inside the mask (a low
  percentile, not the median: between dense particles the ring is mostly
  other particles). `contrast` $=$ peak $-$ background;
* the region is cut at the particle's Voronoi cell (each pixel belongs to
  the nearest centre), so two touching particles cannot merge, and at
  $1.5\, r_\text{max}$ from the centre;
* a particle whose centre is not above the half level gets NaN and is
  dropped.

For a sharp projected sphere of radius $R$ the chord profile is
$\sqrt{R^2 - \rho^2}$; it falls to half its peak at

$$
\rho_{1/2} = \frac{\sqrt 3}{2} R = 0.87\, R ,
$$

so a well-resolved sphere reads 10 to 15 % small. Blur works the other
way: a particle only a few pixels across is sized by the point-spread
function and reads large. The value is reported as measured, without a
model correction, because the calibration ({doc}`calibration`) reproduces
the same bias in the simulated images.

## 5. Measure (`ferrospot.stats`)

`neighbours` finds, for every particle, the nearest other centre with a
k-d tree ($k = 2$):

$$
d_\text{NN}(i) = \min_{j \ne i} |\mathbf x_i - \mathbf x_j|, \qquad
g(i) = d_\text{NN}(i) - r_i - r_j .
$$

$g < 0$ means the two overlap *in projection*; it says nothing about
touching in 3D.

**Border correction.** A particle enters the spacing statistics only if
it is farther from the mask edge (or the image edge, the mask is padded
by one pixel) than from its nearest neighbour, $e_i > d_\text{NN}(i)$;
otherwise a closer neighbour could sit just outside the counted region.
This is the "reduced sample" (minus-sampling) estimator; it is unbiased
but discards the particles near the edge.

`spacing_stats` reports the count $n$, the counted area $A$ (mask pixels
times $s^2$), the areal density $N_A = n / A$, and mean, standard
deviation and 10 / 50 / 90 percentiles of the diameter $2r$, of
$d_\text{NN}$ and of $g$ over the border-corrected particles, the
overlap fraction (share with $g < 0$), and

$$
\text{CE} = \frac{\overline{d_\text{NN}}}{\;0.5 / \sqrt{N_A}\;},
$$

the Clark-Evans ratio {cite:p}`ClarkEvans1954`: 1 for a random (Poisson)
pattern of the same density, above 1 more regular (hard-core exclusion),
below 1 clustered.

With a section thickness $t$ it adds the Poisson 3D estimate. A sphere is
seen when its centre lies within the section or less than one radius
outside it, so

$$
N_A = n_V\,(t + \bar d), \qquad
\overline{d_\text{NN}^{3D}} = \Gamma(4/3)\left(\tfrac{4\pi}{3} n_V\right)^{-1/3}
= 0.554\, n_V^{-1/3}, \qquad
\phi = n_V \left\langle \tfrac{4}{3}\pi r^3 \right\rangle ,
$$

the Hertz–Chandrasekhar mean nearest-neighbour distance of a Poisson
process {cite:p}`Chandrasekhar1943`. It ignores that spheres cannot
overlap and that the detector misses particles, so it is a first guess;
the forward-model calibration replaces it.

`depth_profile` bins the same per-particle numbers by the depth of each
centre below the vacuum (bin width `depth_bin_nm`), with the neighbours
searched among all particles in the mask, not only those in the bin.

## Configs, CLI and dashboard

`ferrospot.config.analyse` runs the five steps from one JSON config per
image (source, scale, mask, detect, expected values, and a free-text
note for every choice); `ferrospot run` and `ferrospot table` run and
tabulate them; the dashboard does the same interactively and writes a
config to start from. The calibration ({doc}`calibration`) reuses
`analyse` and the same `detect` on simulated images.
