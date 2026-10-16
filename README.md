<a id="readme-top"></a>

<br />
<div align="center">
  <h3 align="center">ferrospot</h3>

  <p align="center">
    Spots Fe: iron nanoparticles in TEM / HAADF-STEM images, their sizes and spacings
  </p>
</div>

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li><a href="#getting-started">Getting Started</a></li>
    <li><a href="#running-it-yourself">Running it yourself</a></li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#library">Library</a></li>
        <li><a href="#configs-and-command-line">Configs and command line</a></li>
        <li><a href="#dashboard">Dashboard</a></li>
        <li><a href="#reading-the-numbers">Reading the numbers</a></li>
        <li><a href="#from-2d-to-3d-forward-model-calibration">From 2D to 3D: forward-model calibration</a></li>
      </ul>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



## About The Project

Space-weathered rims of lunar and asteroidal grains hold nanophase iron
(npFe0): bright spheres in HAADF-STEM and dark-field images, dark ones in
bright field, a few nm to a few tens of nm across. `ferrospot`

* finds them with a multi-scale Laplacian-of-Gaussian blob detector
  (mixed sizes) or Difference-of-Gaussians local maxima (dense, small
  particles), in dark- or bright-field images;
* sizes every particle by half-maximum region growing, independent of
  the detector;
* restricts the count to a mask: an intensity window, boxes, a depth band
  below the vacuum interface;
* reports nearest-neighbour centre distances and surface gaps with a
  border correction, the share of particles overlapping in projection,
  the Clark-Evans ratio against a random pattern, statistics in depth
  bins, and a rough 3D estimate for a section of known thickness;
* turns the 2D numbers into true 3D spacings by simulating each image:
  3D packings of iron spheres in a section, projected with a
  Beer-Lambert image model, run through the same pipeline;
* reads figures straight from PDFs (embedded images, vector scale bars);
* keeps one JSON config per analysed image, with the reasons for every
  choice, and runs them from the command line.

A Streamlit dashboard tunes the settings per image and keeps notes.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Getting Started

### Prerequisites

* Python 3.13 or newer
* [uv](https://docs.astral.sh/uv/) (optionally [devenv](https://devenv.sh);
  `direnv allow` picks it up)
* For figures in PDFs only: `pdfimages` and `pdftoppm` from poppler-utils
  (in the devenv shell; elsewhere e.g. `nix shell nixpkgs#poppler-utils`
  or `apt install poppler-utils`)

### Installation

```sh
cd ferrospot
uv sync --extra dashboard   # library, dashboard, pytest and ruff
uv run pytest
```

The library needs NumPy, SciPy, Pillow and scikit-image; Streamlit and
Plotly come only with the `dashboard` extra, polydisc (with numba) and
matplotlib with the `simulate` extra (`uv sync --extra dashboard --extra
simulate`) for the 2D → 3D calibration.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Running it yourself

```sh
uv sync --extra simulate --extra dashboard
# Configs that read figures from PDFs need poppler-utils (pdfimages, pdftoppm).
uv run ferrospot run "configs/*.json" --out results      # detect, measure
uv run ferrospot table results                              # Markdown table
# 2D -> 3D calibration, default grid:
uv run ferrospot calibrate "configs/*.json" --out calibration --workers 8
# with the model upgrades (a new --out: saved grids are reused only when
# made with the same options):
uv run ferrospot calibrate "configs/*.json" --out calibration_v2 --workers 8 \
    --alphas 1,1.1,1.25,1.5,2,2.5,3,4,6,9 --bragg --match nndd
#   --depth exp:30 or --depth peak:20:10 for a rim with the grain surface
#   at one face of the section (applied to every panel given)
uv run ferrospot-dashboard                                  # Streamlit UI
```

`--bragg` acts on bright-field panels only; `--match nndd` fits the 2D
NN-distance distribution and both estimates (mean and distribution) are
written and tabulated. Expect roughly 1.5 to 2 times the default run time
with all options (the alpha family adds ~10 states, `--bragg` doubles the
bright-field grid).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Usage

### Library

```python
import ferrospot

img = ferrospot.load_image("grain.tif")                  # grayscale float
nm_per_px = 0.29                                         # metadata, or:
# nm_per_px = ferrospot.measure_scale_bar(img, 20, region=(0, 450, 200, 512))

# Everything but the dark vacuum; depth = distance to it, in nm.
mask, depth = ferrospot.build_mask(img, nm_per_px, vacuum_below=0.05, smooth_nm=2)

p = ferrospot.detect(img, nm_per_px, detector="dog", radius_nm=(0.8, 3),
                     threshold=0.023, mask=mask, bright_field=False)
print(len(p), 2 * p.r_nm.mean())                         # half-maximum diameter

s = ferrospot.spacing_stats(p, mask=mask, section_thickness_nm=100)
print(s.nn_centre_mean_nm, s.nn_gap_mean_nm, s.overlap_fraction, s.clark_evans)

rows = ferrospot.depth_profile(p, depth, bin_nm=20, mask=mask)
nb = ferrospot.neighbours(p, mask)                       # per particle
p.to_csv("grain_particles.csv", nn_nm=nb["nn_nm"], gap_nm=nb["gap_nm"])
rgb = ferrospot.overlay(img, p, mask)                    # look before you trust
```

* `detect(detector="log")` searches 20 scales between the two radii;
  `"dog"` is one band-pass tuned to the smallest radius, with local
  maxima at least that far apart, and separates dense particles that LoG
  merges. `bright_field=True` inverts the image first. The image is
  rescaled (0.5 to 99.9 percentile inside the mask) so `threshold` means
  roughly the same across images, but it still sets the count: try a few
  values and look at the overlay.
* Sizes: the equivalent-circle radius of the connected region above half
  of the particle's peak over its local background, inside its Voronoi
  cell (so touching neighbours do not merge) and inside the mask.
* `build_mask` combines `roi`, `exclude` boxes, `window` (intensity range,
  holes filled, `erode_nm` from its edges) and a vacuum level
  (`vacuum_below` for dark-field, `vacuum_above` for bright-field) with
  an optional `depth_band_nm`. Intensities are in the loaded image's
  units.
* `ferrospot.pdf.extract_image(pdf, page, index)`,
  `render_page(pdf, page, dpi)` and `image_ppi(pdf, page, index)` reach
  figures in PDFs; image numbers are those of `pdfimages -list -f page
  -l page`.

### Configs and command line

One JSON file per analysed image records the source (an image file, or
PDF + page + embedded image + crop), the scale and how it was measured
(`nm_per_px`, a bar measured by hand, a burned-in bar in a box, or a
vector bar in a box of a page render), the mask and the detection
settings, with free-text `scale_note`, `mask_note` and `threshold_note`
for every choice that is not obvious; see the docstring of
`ferrospot.config`.

```sh
uv run ferrospot run "configs/*.json" --out results  # per config: summary JSON,
                                                     # particles CSV, overlay PNG
uv run ferrospot table results                       # Markdown summary table
```

The dashboard's *Config JSON* button writes a starting point.

### Dashboard

```sh
uv run ferrospot-dashboard       # or: streamlit run src/ferrospot/dashboard.py
```

Pick an image from a folder or upload one, set nm/px (or measure a bar),
choose LoG or DoG, dark or bright field, radius range and threshold, and
build the mask (ROI, excluded boxes, intensity window with erosion,
vacuum level with depth band). The page shows the half-maximum circles
(red markers count in the spacing statistics, orange ones fail the
border correction, cyan is the mask edge; hover for size, NN distance
and gap), the diameter, NN distance, gap, overlap and Clark-Evans ratio,
all statistics, and a depth profile when a vacuum level is set.
Download the particles as CSV or the settings as a config. The comment
box appends a timestamped entry with the current parameters and
statistics to `<notes folder>/<image-stem>_<YYYY-MM-DD>.md` (default
folder `./ferrospot-notes/`); earlier notes for the image are shown
below it.

### Reading the numbers

The assumptions are also in the docstring of `ferrospot.stats`.

* The image is a projection of a section of thickness *t*: particles at
  different depths overlap, so 2D nearest-neighbour distances are shorter
  than the 3D ones and 2D gaps can be negative. "Overlapping" means
  overlapping in projection, not touching.
* Detection pushes the other way: merged neighbours and missed small or
  faint particles make particles look larger and farther apart.
* Half-maximum diameters of a sharp projected sphere are 0.87 of the
  true diameter; particles only a few pixels across are sized by the
  blur and read large.
* Border correction: a particle enters the spacing statistics only if it
  is farther from the mask or image edge than from its nearest neighbour.
* Clark-Evans ratio: mean NN distance over 0.5/sqrt(density), its value
  for a random (Poisson) pattern; above 1 is more regular (e.g.
  exclusion), below 1 clustered.
* 3D estimate, random (Poisson) placement: a particle shows up when its
  centre lies within the section or less than a radius outside it, so
  n<sub>V</sub> = N<sub>A</sub> / (*t* + *d*<sub>mean</sub>); the mean 3D
  nearest-neighbour distance is 0.554 n<sub>V</sub><sup>-1/3</sup>, the
  gap that minus *d*<sub>mean</sub>, and the volume fraction
  n<sub>V</sub> ⟨4/3 π r³⟩. Poisson ignores that particles cannot overlap,
  so treat it as a rough guide once the fill reaches a few percent. It
  also assumes every particle is found; the forward-model calibration
  below does not.

### From 2D to 3D: forward-model calibration

The image of a 50 to 100 nm section is a projection: particles at
different depths overlap, so a mean 2D NN distance of 3 nm can come from
particles that are several radii apart in 3D. `ferrospot.calibrate`
simulates each panel and reads the true 3D state off its measured 2D NN
distance (method and assumptions in
[docs/calibration.md](docs/calibration.md)):

```python
import ferrospot
from ferrospot import calibrate, simulate

# 3D packing -> image, by hand:
sec = simulate.section(87, 100, simulate.radius_spec(1.2, 1.25), gap_nm=5, seed=0)
img = simulate.simulate_image(sec, 0.29, mode="haadf", psf_fwhm_nm=0.6, noise=0.03)
print(simulate.truth(sec))          # true 3D NN gap, fill, number density

# A panel: simulate its imaging conditions over a grid of 3D states
cfg = ferrospot.load_config("configs/my_panel.json")
res = ferrospot.analyse(cfg, "configs")
cal = calibrate.calibrate_panel(cfg, res, workers=2)
est = calibrate.estimate_3d(res["stats"], cal)
print(est["gap3d_r"], est["gap3d_r_lo"], est["gap3d_r_hi"], est["completeness"])
```

```sh
uv run ferrospot calibrate "configs/*.json" --out calibration --workers 2
```

writes per panel the calibration (reused next time), the estimate, a
plot of the calibration curves and a real-versus-simulated comparison,
and a Markdown table. Optional model upgrades, all off by default: a
radius-scaled hard core (`--alphas`, `simulate.section(...,
packing="alpha", alpha=1.5)`), random Bragg visibility in bright field
(`--bragg`, `simulate_image(..., bragg=True)`), a depth profile of the
number density (`--depth`, `section(..., depth={"kind": "exp",
"scale_nm": 30})`) and matching the whole NN-distance distribution
(`--match nndd`, `estimate_3d(..., match="nndd")`); see
[docs/calibration.md](docs/calibration.md#model-options). The dashboard's *Simulate a section* page does the
same for one section interactively: set thickness, sizes, 3D gap,
imaging and detection, and compare what is measured with what is there.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Roadmap

- [x] Forward-model calibration of the 2D → 3D bias
      ([docs/calibration.md](docs/calibration.md)).
- [x] Radius-scaled hard core, Bragg visibility, depth profiles and
      NN-distribution matching in the calibration (options).
- [ ] 2D pair-correlation g(r) as a further matching criterion.
- [ ] Depth-graded sizes in the simulated rims, and calibrating
      depth bins rather than whole panels.
- [ ] Use the same synthetic images to measure the size and count bias of
      each detector and threshold (the threshold has no plateau in most
      real panels).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



## License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for more
information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Contact

Mirza Arnaut - mirza.arnaut@tu-dortmund.de

<p align="right">(<a href="#readme-top">back to top</a>)</p>
