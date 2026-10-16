# Getting started

## Install

```sh
pip install ferrospot                      # library: numpy, scipy, pillow, scikit-image
pip install "ferrospot[dashboard]"         # + streamlit, plotly
pip install "ferrospot[simulate]"          # + polydisc (numba), matplotlib: 2D → 3D calibration
```

From source, with [uv](https://docs.astral.sh/uv/):

```sh
git clone https://github.com/arunoruto/ferrospot && cd ferrospot
uv sync --extra dashboard --extra simulate   # also installs pytest, ruff and sphinx
uv run pytest -m "not data"
uv run sphinx-build -b html docs docs/_build/html   # this site
```

With [devenv](https://devenv.sh) (`direnv allow`), `poppler-utils` is on
the path for figures in PDFs, and `devenv tasks run docs:build` builds
the site.

## Library

```python
import ferrospot

img = ferrospot.load_image("grain.tif")                   # 2D float array
s = 0.29                                                  # nm per pixel (metadata), or
# s = ferrospot.measure_scale_bar(img, 20, region=(0, 450, 200, 512))

mask, depth = ferrospot.build_mask(img, s, vacuum_below=0.05, smooth_nm=2)
p = ferrospot.detect(img, s, detector="dog", radius_nm=(0.8, 3), threshold=0.023, mask=mask)
st = ferrospot.spacing_stats(p, mask=mask, section_thickness_nm=100)
print(len(p), st.d_mean_nm, st.nn_centre_mean_nm, st.nn_gap_mean_nm, st.clark_evans)
```

{doc}`pipeline` explains each step and its parameters.

## Configs and command line

One JSON file per image records the source, scale, mask, detection
settings, and a note for every choice (the docstring of
`ferrospot.config` lists the keys):

```sh
ferrospot run "configs/*.json" --out results   # summary JSON, particles CSV, overlay PNG
ferrospot table results                        # Markdown table
ferrospot calibrate "configs/*.json" --out calibration --workers 4   # 2D → 3D, see Calibration
```

## Dashboard

```sh
ferrospot-dashboard            # or: streamlit run src/ferrospot/dashboard.py
```

Pick or upload an image, set the pixel size, tune the mask and the
detector while watching the overlay, download the particles as CSV or the
settings as a config, and keep timestamped notes per image. A second page
simulates a section and shows what the pipeline measures in it.
