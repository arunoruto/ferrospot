# Dashboard: current state and plan

The dashboard (`ferrospot.dashboard`, Streamlit) has two pages chosen by a
sidebar radio: *Analyse an image* and *Simulate a section*. Everything
below is the plan for the next version; nothing here is implemented yet.
Facts about Streamlit and Plotly were checked against the installed
versions (Streamlit 1.66, Plotly 7.1) and their documentation on
2026-10-16.

## What stays

* **Streamlit.** The alternatives (Bokeh/Panel, Dash, napari) would be a
  rewrite for interactions that Streamlit now gives for free.
* **Plotly for the image.** `st.plotly_chart(fig, on_select="rerun")`
  returns clicked points, box selections and lasso selections to Python
  since Streamlit 1.35, so "click a particle and something custom
  happens" needs no extra component. `px.imshow` axes are in pixel
  units, so every returned coordinate is an exact pixel coordinate
  whatever the zoom.
* `@st.cache_data` on the analysis with the image array as an argument:
  hashing a 2000 × 2000 float image costs about 25 ms (shape, dtype and
  a 100 000-element sample), so no custom hash is needed.

## Layout: one entry point, pages as functions

Use `st.navigation` + `st.Page` with callables rather than a `pages/`
directory. A `pages/` folder must sit next to the entry-point file and is
found by filesystem lookup, which is brittle for a package installed in
`site-packages` and launched through the `ferrospot-dashboard` console
script. Callables need no files on disk, order and titles are explicit,
and the docs recommend them "for maximum flexibility".

```
src/ferrospot/dashboard/
    __init__.py      # main(): streamlit run on app.py (console script, unchanged)
    app.py           # set_page_config, sidebar image picker, st.navigation([...]).run()
    analyse.py       # def page(): ...
    simulate.py      # def page(): ...
    label.py         # def page(): ...
    common.py        # _load, _figure, box parsing, session-state keys
```

```python
# app.py
st.set_page_config(page_title="ferrospot", layout="wide")
name, image = pick_image(st.sidebar)          # shared by every page
st.navigation([
    st.Page(analyse.page,  title="Analyse",  icon=":material/analytics:", url_path="analyse", default=True),
    st.Page(simulate.page, title="Simulate", icon=":material/science:",   url_path="simulate"),
    st.Page(label.page,    title="Label",    icon=":material/edit:",      url_path="label"),
]).run()
```

The entry point reruns on every interaction, pages only when selected;
widgets shared by all pages (image picker, nm/px) live in the entry point
("Option 1 (preferred)" in the Streamlit docs on widgets across pages).
`st.session_state` is one dict per browser session and is shared by all
pages, so the image, pixel size, boxes and manual edits can be set on one
page and used on another.

## The four interactions

All four come from one figure built by `_figure` (image, mask outline,
particle markers) and one `st.plotly_chart(..., on_select="rerun",
key=...)` call. What `on_select` returns, verified in the Streamlit
frontend bundle:

* `selection.points`: list of `{curve_number, point_number, point_index,
  x, y, customdata, ...}`. A single click on a marker is a selection
  (Streamlit sets `clickmode = "event+select"` when `"points"` is in
  `selection_mode`); clicking it again deselects (empty selection).
* `selection.box`: list of `{x: [x0, x1], y: [y0, y1]}` from box select;
  Shift-drag adds further boxes.
* `selection.lasso`: list of `{x: [...], y: [...]}`.
* **Not** returned: shapes drawn with Plotly's `drawrect` modebar
  buttons. They stay in the browser (the frontend subscribes to select,
  click and deselect, not to relayout). Do not plan on them.

### 1. ROI and exclude boxes: Plotly box select

```python
fig.update_layout(dragmode="select")
ev = st.plotly_chart(fig, key="roi", on_select="rerun", selection_mode=("box",))
for b in ev.selection.box:
    x0, x1 = sorted(b["x"]); y0, y1 = sorted(b["y"])
    st.session_state.boxes.append((int(x0), int(y0), int(x1), int(y1)))
```

A radio chooses whether the next box is the ROI or an exclusion. Boxes
are kept in session state, drawn back with `fig.add_shape(type="rect")`
(display only) and listed with a delete button each; they feed
`build_mask(roi=..., exclude=...)` exactly as the typed boxes do now, and
go into the downloaded config.

### 2. Click to add or remove particles

A click on empty image area produces no event: image traces do not emit
point selections. The cheap trick is an invisible grid trace so every
click hits some point:

```python
step = 8  # px
gx, gy = np.meshgrid(np.arange(0, w, step), np.arange(0, h, step))
fig.add_scatter(x=gx.ravel(), y=gy.ravel(), mode="markers",
                marker={"size": step, "opacity": 0}, hoverinfo="skip", name="grid")
fig.add_scatter(x=p.x_nm / s, y=p.y_nm / s, mode="markers", name="particles")
ev = st.plotly_chart(fig, key="label", on_select="rerun", selection_mode=("points",))
for pt in ev.selection.points:
    if pt["curve_number"] == PARTICLES:
        edits.append(("remove", pt["point_index"]))
    else:
        edits.append(("add", pt["x"], pt["y"]))
```

The grid is about 62 000 points for a 2000 × 2000 image at 8 px, which
Plotly handles; the click resolution is the grid step, good enough to
place a particle whose size is then measured by `half_max_radii` from
the image (so a manually added particle gets the same half-maximum
radius as a detected one). Edits are a list of tuples in session state;
"Undo" pops the last one; the effective particle set is derived from the
detections plus the edits on every run (pure and cheap). Edits are saved
in the config as `"manual": {"add": [[x, y], ...], "remove": [[x, y], ...]}`
in pixels, and `analyse` applies them after detection, so a CLI run
reproduces the dashboard result. Deduplicate events: Streamlit resends an
empty selection on deselect and double-click, so store the last processed
event and ignore repeats.

### 3. Two clicks to measure the scale bar

Same mechanism: accumulate two clicks on the grid trace, then
`nm_per_px = bar_nm / dist(p1, p2)`. The zero-trick alternative is
`streamlit-image-coordinates` (MIT, maintained, Streamlit ≥ 1.34) with
`click_and_drag=True`, which returns `{x1, y1, x2, y2}` for one drag; it
shows the image as a plain `<img>`, so for the scale bar alone it is the
least code. Keep the existing threshold-based helper as a third option.

### 4. Click a particle for its profile and neighbours

`pt["point_index"]` on the particle trace indexes into `Particles`. Show
next to the chart: diameter, contrast, NN distance and gap (already in
`customdata`), the intensity profile through the centre along x and y
with the half-maximum level and the fitted outline, and lines to the
nearest neighbours drawn as a short-lived shape on the figure.

## Label page

A page for building ground truth by hand, so the detectors can be scored
on real images, not only on simulations:

* starts from the current detections (or none), with the add/remove
  mechanism above and the same radius measurement;
* a size slider and an "accept size" toggle per particle for the cases
  where half-maximum sizing fails (touching particles, bright-field
  diffraction contrast);
* saves `<image-stem>_labels.csv` (`x_nm, y_nm, r_nm, source`) next to
  the notes, and loads it back;
* a score panel: precision, recall and size error of the current detector
  settings against the labels (match within one radius, as
  `calibrate._match` does), so the threshold can be tuned against the
  labels instead of by eye.

## Keeping it responsive

* The figure JSON carries the whole image on every rerun: 5.7 MB for a
  2000 × 2000 8-bit grey noise image, 16 MB for RGB, 54 MB for a float
  heatmap. Send the overlay as `uint8` and let real images compress;
  there is no incremental update API in `st.plotly_chart`.
* Wrap chart and click handling in `@st.fragment`, so a click reruns the
  fragment, not the page; cache the base figure dict with
  `@st.cache_data` and only append the cheap marker traces per run.
* `theme=None` skips Streamlit's rewrite of the figure JSON.

## Components considered and not chosen

| Package | Verdict |
|---|---|
| `streamlit-drawable-canvas` 0.13 (MIT, needs Streamlit ≥ 1.53) | Works again since 0.10; gives free-drawn, editable, labelled rectangles and returns boxes in image pixels. Only worth it if many persistent labelled boxes are needed; the particle overlay would have to be baked into its background image. |
| `streamlit-image-coordinates` 0.4 (MIT) | Good for the scale bar; one click or one drag per event, no overlay. |
| `streamlit-image-annotation` 0.8 (Apache-2.0) | Bounding boxes and points from a file path; a second widget beside the Plotly figure. |
| `streamlit-plotly-events` | Obsolete, broken on Streamlit ≥ 1.42; replaced by `on_select`. |
| `streamlit-bokeh` (official) | Returns nothing to Python: tap and box-select events do not reach the app. |
| `streamlit-label-kit` | GPL-2.0-only. |
| `streamlit-cropper`, `-cropperjs`, `-img-label`, `-image-zoom` | Stale or single-box only. |
