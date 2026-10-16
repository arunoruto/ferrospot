# From 2D to 3D: forward-model calibration

Every spacing ferrospot measures is a 2D projection through a section
50 to 100 nm thick. Particles at different depths overlap in projection,
so the measured nearest-neighbour (NN) distance, surface gap and overlap
fraction understate the true 3D spacing; merged and missed particles push
the other way. What matters for optical coupling between iron spheres is
the true 3D surface gap in units of the radius (pair coupling roughly
34 % at 0.1 r, 14 % at 0.25 r, 5 % at 0.5 r, 1.5 % at 1 r, 0.2 % at 2 r).

`ferrospot.calibrate` does not correct the 2D numbers analytically. It
simulates sections with known 3D packings under the imaging conditions of
one panel, runs exactly the same detection and statistics on them, and
reads off which 3D states give the panel's measured 2D NN distance.

## Model

**Packings** (`simulate.section`). [polydisc](https://pypi.org/project/polydisc/)
fills a box, periodic in x and y, with spheres with lognormal radii by
random sequential addition to saturation, with a minimum surface gap.
The box is thicker than the section by a margin on both sides, so its
walls do not order the spheres inside; spheres reaching out of the
section are cut by it (as in a FIB lamella). Two families of states:

* *regular*: saturated at a minimum gap of 0.25 to 16 median radii. The
  gap fixes the density; the spacing is more even than random.
* *random*: a saturated packing without a gap, thinned at random to the
  densities of the regular states and below. (polydisc's `max_points`
  cannot be used for sparse states: it stops a growth that starts from
  one point and leaves a compact cluster.)

At the same density the random family has a broader NN distribution and
a smaller mean 3D gap; the difference between the two is part of the
reported uncertainty.

**Model options** (off by default; the default grid and results below
are unchanged). <a id="model-options"></a>

* *Radius-scaled hard core* (`section(..., packing="alpha", alpha=a)`,
  CLI `--alphas 1,1.25,1.5,2,3`): centres of spheres i and j at least
  `a (r_i + r_j)` apart (`a = 1`: touching allowed), lognormal radii,
  saturated. Implemented as a polydisc packing of spheres `a` times
  larger shrunk back, so periodic in x and y like the others. With
  `--alphas` the grid gets a family `alpha` (saturated at each `a`,
  dotted in the plots); the estimate then also reports `alpha` with its
  range (from that family; the gap and fill ranges span all families).
  A fixed gap in nm spaces small particles relatively further apart than
  large ones; a radius-scaled one does not.
* *Bragg visibility* (`simulate_image(..., bragg=True | dict)`, CLI
  `--bragg`): per sphere a random factor on its `k_Fe`: `k_strong = 4`
  for a share `p_strong = 0.25` (oriented for strong diffraction),
  `k_mid = 2` for `p_mid = 0.25`, `k_weak = 1` (mass-thickness only) for
  the rest (`simulate.BRAGG`; a dict overrides any of them). Seeded with
  the image seed, independently of the noise. Bright field only: in HAADF
  the excess over 1 is multiplied by `haadf` (default 0; e.g. 0.1 for a
  small residual diffraction effect at low camera length). In the
  calibration (bright-field panels only) every state is simulated with
  and without it (families `random+bragg`, `regular+bragg`, ...), a
  nuisance folded into the range; the pilot (noise, contrast, size
  matching) runs without it.
* *Depth profile* (`section(..., depth=...)`, CLI `--depth exp:30` or
  `--depth peak:20:10`): number density along z for a section whose face
  `z = 0` is the grain surface. Spheres are kept with probability
  `exp(-z / L)` or `exp(-(z - z0)^2 / 2 w^2)`, none above the surface,
  after packing (so local packings stay hard-core; the density falls by
  thinning). Applies to every state and the pilot.
* *Distribution matching* (`estimate_3d(..., match="nndd")`, CLI
  `--match nndd`): every simulated image records 19 quantiles (5 to
  95 %) of its border-corrected 2D NN distances, the panel too. The
  distance of two distributions (`calibrate.nndd`) is the mean absolute
  difference of the quantiles over the reference mean (a Wasserstein-1
  distance over the central 90 %, relative: 0.1 for a 10 % scale
  change). Per curve the estimate is the state of least distance;
  consistent are those within 10 % (plus the seed scatter) of that
  least distance. The least distance itself is reported (`nndd`, the
  "misfit" in the table): well above 0.1 means no simulated state has
  the panel's NN distribution, whatever its mean. The CLI writes both
  estimates (`other` in `<name>_estimate.json`) and the table adds a
  "from NNDD" column. Not implemented: the pair-correlation g(r).

Caveats: `--bragg` changes the particle contrast after the pilot matched
it, so the `+bragg` families have darker particles than the panel's
median; the alpha family and `--depth` change the densest reachable
state (alpha 1 saturates at the same fill as the `regular` gap 0, a
surface-peaked profile lowers the mean density); a saved calibration is
reused only if made with the same options, so use a new `--out` to keep
the default results.

**Image** (`simulate.simulate_image`), projected-thickness Beer-Lambert.
Per pixel the electron path crosses `L` nm of iron (the sum of the
sphere chords inside the section, sampled on a 3x finer grid) and
`t - L` nm of matrix; with `tau = k_m (t - L) + k_Fe L`,

* bright field: `I = exp(-tau)` (scattering out of the objective
  aperture; particles dark);
* HAADF / dark field: `S = 1 - exp(-tau)` (scattering onto the annular
  detector; particles bright, saturating for thick heavy columns).

At 100 nm both are nearly linear in mass-thickness, and since the image
is stretched afterwards only `k_Fe / k_m` matters: 4 for HAADF
(incoherent scattering grows as Z^1.7; Fe, Z = 26, against a silicate of
mean Z about 11 at a similar atom density) and 3 for bright field
(mass-thickness). The matrix thickness varies smoothly by a few percent
(a Gaussian random field over an eighth of the image): real sections are
uneven, and that variation, not the particles, sets much of a published
figure's grey range. Then a Gaussian blur (FWHM default 2 px: probe and
figure resampling), averaging into the panel's pixels, a percentile
stretch (0.5 to 99.9 %), white Gaussian noise and 8-bit quantisation.

**Matching a panel** (`calibrate.calibrate_panel`). From the real panel:
pixel size, mode, detector, radius range and threshold (its config), the
mean and spread of the half-maximum diameters, the pixel-scale noise (the
image minus its 1-px Gaussian smoothing, robust standard deviation, after
the same normalisation as detection) and the median detected particle
contrast. The radii are lognormal with the measured spread, cut off below
the detector's smallest radius (`radius_nm[0]`): the calibrated
population is what the panel's settings count, and its 3D spacing is
that of those particles. (Without the cut a broad measured spread puts
thousands of invisible sub-nm spheres into the model, and the "3D gap"
becomes theirs.) A pilot loop on a random section adjusts, in turn, the noise
level, the thickness variation (contrast), the median radius (diameter)
and the density (NN distance) until the simulated values match. Then the
grid: both families, 10 regular gaps plus 2 sparser random states, three
section thicknesses, three seeds: about 200 images per panel.

**Inverse** (`calibrate.estimate_3d`). Per thickness and family the seed
means form a curve of measured 2D NN distance against true 3D gap. States
whose simulated NN distance lies within 10 % of the panel's (the change a
small threshold change makes), plus
the seed scatter, are consistent. The estimate is the mean of the two
families at the middle thickness; the range spans all consistent states
of all thicknesses and both families. Where the curve is flat (detection
saturated: the detector cannot separate particles closer than about
`r_min`) the range runs to the end of the grid. The simulated count
density, overlap fraction and Clark-Evans ratio at the estimate are not
fitted; their ratio to the measured values is a check of the model.

True 3D quantities, over the spheres with their centre in the section and
neighbours searched in the whole slab: mean NN centre distance, mean NN
surface gap (in nm, and over the mean radius), the share of spheres whose
gap is below their own radius ("gap < r"), the volume fill and number
density. *Completeness* is the share of spheres reaching into the section
that a detection sits on (within its radius or 1 px).

## Assumptions

* Kinematic, projected-thickness contrast. Diffraction contrast of
  crystalline iron (in bright field some particles are much darker than
  others for that reason) only as the optional random Bragg factor,
  not from orientations; no phase or Fresnel contrast, no dynamical
  scattering or channelling.
* Particles smaller than the detector's smallest radius are not
  modelled; if a population of them exists, the true spacing of all
  iron is smaller than reported.
* Spheres only, lognormal sizes, uniform in depth (unless `--depth`)
  and laterally: no depth-graded sizes, no clustering. Real clustering would put more particles within
  one radius than either family.
* Uniform section thickness apart from the few-percent variation; the
  thickness range is a prior set per config (`thicknesses_nm`): the value
  from the sample-preparation methods +-20 % where it is stated, else
  50, 75 and 100 nm.
* The whole panel is one population: a calibration per panel, not per
  depth bin.
* PSF 2 px FWHM for every panel; the size matching absorbs part of an
  error in it.

## Usage

```sh
uv sync --extra simulate
uv run ferrospot calibrate "configs/*.json" --out calibration --workers 2
# with the options (see Model options):
uv run ferrospot calibrate "configs/*.json" --out calibration_v2 --workers 8 \
    --alphas 1,1.1,1.25,1.5,2,2.5,3,4,6,9 --bragg --match nndd
```

writes `<name>_calibration.json` (the simulated grid; reused next time),
`<name>_estimate.json`, `<name>_calibration.png` (the curves),
`<name>_compare.png` (a crop of the real panel next to a simulation at the
inferred state; keep it local, the panels are copyrighted) and
`table.md`. About 2 to 10 minutes per panel on two cores.

```python
import ferrospot
from ferrospot import calibrate, simulate

cfg = ferrospot.load_config("configs/my_panel.json")
res = ferrospot.analyse(cfg, "configs")
cal = calibrate.calibrate_panel(cfg, res, workers=2)   # or thicknesses_nm=[...]
est = calibrate.estimate_3d(res["stats"], cal)
# options: calibrate_panel(..., alphas=(1, 1.5, 2, 3), bragg=True,
#                          depth={"kind": "exp", "scale_nm": 30})
# estimate_3d(res["stats"], cal, match="nndd")   # distribution, not mean
est["gap3d_r"], est["gap3d_r_lo"], est["gap3d_r_hi"]
calibrate.plot_calibration(cal, est).savefig("curves.png")
calibrate.plot_strip(cal).savefig("strip.png")

# Or one section by hand
sec = simulate.section(87, 100, simulate.radius_spec(1.2, 1.25), gap_nm=5, seed=0)
img = simulate.simulate_image(sec, 0.29, mode="haadf", noise=0.03)
simulate.truth(sec)
sec = simulate.section(87, 100, simulate.radius_spec(1.2, 1.25), packing="alpha",
                       alpha=1.5, depth={"kind": "peak", "z0_nm": 20, "w_nm": 10})
img = simulate.simulate_image(sec, 0.29, mode="bf", bragg={"p_strong": 0.3})
```

The dashboard's *Simulate a section* page shows one simulated section
with its detections, the true 3D and the measured 2D numbers side by
side.
