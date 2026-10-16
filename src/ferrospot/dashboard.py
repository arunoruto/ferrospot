"""Streamlit dashboard: pick an image, tune mask and detection, save notes;
or simulate a section and see what the pipeline measures in it.

Run with ``ferrospot-dashboard`` (installed with the ``dashboard`` extra) or
``streamlit run src/ferrospot/dashboard.py``.
"""

import io
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import plotly.express as px
import streamlit as st

import ferrospot

SUFFIXES = {".png", ".tif", ".tiff", ".jpg", ".jpeg"}


def main():
    """Entry point of ``ferrospot-dashboard``: run this file with Streamlit."""
    from streamlit.web import cli

    sys.argv = ["streamlit", "run", __file__, *sys.argv[1:]]
    sys.exit(cli.main())


@st.cache_data
def _load(path_or_bytes):
    if isinstance(path_or_bytes, bytes):
        path_or_bytes = io.BytesIO(path_or_bytes)
    return ferrospot.load_image(path_or_bytes)


@st.cache_data
def _analyse(image, nm_per_px, mask_kw, detect_kw, thickness, bin_nm):
    mask, depth = ferrospot.build_mask(image, nm_per_px, **mask_kw)
    p = ferrospot.detect(image, nm_per_px, mask=mask, **detect_kw)
    stats = ferrospot.spacing_stats(p, mask=mask, section_thickness_nm=thickness)
    profile = [] if depth is None else ferrospot.depth_profile(p, depth, bin_nm, mask)
    return mask, p, ferrospot.neighbours(p, mask), stats, profile


def _box(text):
    """``"x0,y0,x1,y1"`` to a tuple of ints; empty text to None."""
    return tuple(int(v) for v in text.split(",")) if text.strip() else None


def _measure(image, bar_nm, box, dark):
    st.session_state.nm_per_px = ferrospot.measure_scale_bar(
        image, bar_nm, _box(box), dark
    )
    st.session_state.bar_px = bar_nm / st.session_state.nm_per_px


def _figure(image, p, mask, nb):
    fig = px.imshow(ferrospot.overlay(image, p, mask))
    s = p.nm_per_px
    fig.add_scatter(
        x=p.x_nm / s,
        y=p.y_nm / s,
        mode="markers",
        marker={"size": 4, "color": np.where(nb["counted"], "red", "orange")},
        customdata=np.column_stack([2 * p.r_nm, nb["nn_nm"], nb["gap_nm"]]),
        hovertemplate="d = %{customdata[0]:.2f} nm<br>NN = %{customdata[1]:.2f} nm"
        "<br>gap = %{customdata[2]:.2f} nm<extra></extra>",
    )
    fig.update_layout(height=750, showlegend=False, margin={"l": 0, "r": 0, "t": 0})
    return fig


def _source(sb):
    """Sidebar image picker; returns (name, image) or None."""
    if sb.radio("Image source", ["Folder", "Upload"], horizontal=True) == "Folder":
        folder = Path(sb.text_input("Folder", "."))
        files = sorted(f for f in folder.glob("*") if f.suffix.lower() in SUFFIXES)
        if not files:
            st.info(f"No PNG/TIFF/JPEG images in {folder.resolve()}.")
            return None
        chosen = sb.selectbox("File", files, format_func=lambda f: f.name)
        return chosen.name, _load(str(chosen))
    upload = sb.file_uploader("Image", type=[s[1:] for s in SUFFIXES])
    if upload is None:
        st.info("Upload an image.")
        return None
    return upload.name, _load(upload.getvalue())


def _mask_controls(sb):
    """Sidebar mask settings, as keyword arguments of build_mask."""
    sb.subheader("Mask")
    kw = {"roi": _box(sb.text_input("ROI x0,y0,x1,y1 (px; empty = all)"))}
    exclude = sb.text_input("Exclude boxes (x0,y0,x1,y1; x0,y0,x1,y1; ...)")
    kw["exclude"] = [_box(b) for b in exclude.split(";") if b.strip()]
    kw["smooth_nm"] = sb.number_input("Smoothing for thresholds (nm)", 0.0, value=1.0)
    if sb.checkbox("Intensity window"):
        kw["window"] = sb.slider("Window (image units)", 0.0, 1.0, (0.1, 0.9), 0.01)
        kw["erode_nm"] = sb.number_input("Erode window edges (nm)", 0.0, value=0.0)
    vacuum = sb.radio("Vacuum", ["none", "dark (below)", "bright (above)"])
    if vacuum != "none":
        level = sb.slider("Vacuum level (image units)", 0.0, 1.0, 0.1, 0.01)
        kw["vacuum_below" if vacuum.startswith("dark") else "vacuum_above"] = level
        if sb.checkbox("Depth band"):
            d0 = sb.number_input("Depth from (nm)", 0.0, value=0.0)
            d1 = sb.number_input("Depth to (nm)", 0.0, value=100.0)
            kw["depth_band_nm"] = (d0, d1)
    return kw


@st.cache_data
def _simulate(size, t, r, gstd, gap_r, keep, s, mode, psf, noise, tsd, seed, kw):
    from ferrospot import calibrate, simulate

    sec = simulate.section(
        size, t, simulate.radius_spec(r, gstd), gap_nm=gap_r * r, keep=keep, seed=seed
    )
    img = simulate.simulate_image(
        sec, s, mode=mode, psf_fwhm_nm=psf, noise=noise, thickness_sd=tsd, seed=seed
    )
    mask = np.ones(img.shape, bool)
    p = ferrospot.detect(img, s, mask=mask, **kw)
    found, _ = calibrate._match(sec, p, img.shape[0] * s)
    stats = ferrospot.spacing_stats(p, mask=mask)
    return (
        img,
        mask,
        p,
        ferrospot.neighbours(p, mask),
        stats,
        simulate.truth(sec),
        found,
    )


def simulate_page(sb):
    """Simulate a section (needs the ``simulate`` extra) and detect in it."""
    sb.subheader("Section")
    t = sb.slider("Thickness (nm)", 5.0, 200.0, 100.0, 5.0)
    r = sb.number_input("Median radius (nm)", 0.2, 100.0, 1.2, 0.1)
    gstd = sb.slider("Geometric std of radii (1 = equal)", 1.0, 2.0, 1.25, 0.05)
    gap_r = sb.select_slider(
        "Minimum 3D gap (radii; saturated packing)",
        [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0],
        4.0,
    )
    keep = sb.slider("Keep (random thinning)", 0.01, 1.0, 1.0, 0.01)
    sb.subheader("Imaging")
    mode = sb.radio("Mode", ["haadf", "bf"], horizontal=True)
    s = sb.number_input("nm per pixel", 0.01, 20.0, 0.29, 0.01, format="%.3f")
    px = sb.slider("Image size (px)", 100, 600, 300, 50)
    psf = sb.number_input("PSF FWHM (nm)", 0.0, 50.0, 2 * s, 0.05)
    noise = sb.slider("Noise (of the grey range)", 0.0, 0.2, 0.03, 0.005)
    tsd = sb.slider("Thickness variation (relative)", 0.0, 0.3, 0.06, 0.01)
    seed = sb.number_input("Seed", 0, 10_000, 0)
    sb.subheader("Detection")
    detector = sb.radio("Detector", ["dog", "log"], horizontal=True)
    kw = {
        "detector": detector,
        "bright_field": mode == "bf",
        "radius_nm": (
            sb.number_input("Min radius (nm)", min_value=0.01, value=0.8),
            sb.number_input("Max radius (nm)", min_value=0.01, value=3.0),
        ),
        "threshold": sb.slider("Threshold", 0.005, 0.5, 0.023, step=0.001),
    }
    try:
        img, mask, p, nb, stats, true, found = _simulate(
            px * s, t, r, gstd, gap_r, keep, s, mode, psf, noise, tsd, seed, kw
        )
    except (ImportError, ValueError) as err:
        st.error(
            f"{err} (install the `simulate` extra; large gaps need a larger image)"
        )
        return
    st.caption(
        "A periodic slab packed with polydisc, projected with Beer-Lambert "
        "absorption (bright field) or scattering (HAADF), blurred, binned, "
        "noised and quantised; then the same detection as for real images."
    )
    cols = st.columns(6)
    cols[0].metric("3D gap / r", f"{true['gap3d_r']:.2f}")
    cols[1].metric("3D gap (nm)", f"{true['gap3d_nm']:.2f}")
    cols[2].metric("3D NN (nm)", f"{true['nn3d_nm']:.2f}")
    cols[3].metric("Gap < r", f"{100 * true['within_1r']:.0f} %")
    cols[4].metric("Fill", f"{100 * true['fill']:.2f} %")
    cols[5].metric(
        "Detected", f"{100 * found:.0f} %", help="of the spheres in the section"
    )
    cols = st.columns(6)
    cols[0].metric("Particles found", f"{stats.count} / {true['n']}")
    for col, label, val in [
        (cols[1], "2D diameter (nm)", stats.d_mean_nm),
        (cols[2], "2D NN (nm)", stats.nn_centre_mean_nm),
        (cols[3], "2D gap (nm)", stats.nn_gap_mean_nm),
        (cols[4], "Overlap in projection", stats.overlap_fraction),
        (cols[5], "Clark-Evans", stats.clark_evans),
    ]:
        col.metric(label, "–" if val is None else f"{val:.2f}")
    st.plotly_chart(_figure(img, p, mask, nb))


def app():
    """The dashboard page."""
    st.set_page_config(page_title="ferrospot", layout="wide")
    sb = st.sidebar
    if (
        sb.radio("Page", ["Analyse an image", "Simulate a section"])
        != "Analyse an image"
    ):
        simulate_page(sb)
        return
    picked = _source(sb)
    if picked is None:
        return
    name, image = picked

    st.session_state.setdefault("nm_per_px", 1.0)
    sb.number_input("nm per pixel", min_value=1e-4, format="%.4f", key="nm_per_px")
    with sb.expander("Scale-bar helper"):
        bar_nm = st.number_input("Bar length (nm)", min_value=0.001, value=100.0)
        box = st.text_input("Search box x0,y0,x1,y1 (px; empty = whole image)")
        dark = st.checkbox("Dark bar on bright background")
        st.button("Measure", on_click=_measure, args=(image, bar_nm, box, dark))
        if "bar_px" in st.session_state:
            st.caption(f"Bar found: {st.session_state.bar_px:.0f} px. Check it.")
    nm_per_px = st.session_state.nm_per_px

    sb.subheader("Detection")
    detector = sb.radio(
        "Detector", ["log", "dog"], horizontal=True, help="LoG: mixed sizes. "
        "DoG: dense small particles, tuned to the minimum radius.",
    )  # fmt: skip
    detect_kw = {
        "detector": detector,
        "bright_field": sb.checkbox("Bright field (dark particles)"),
        "radius_nm": (
            sb.number_input("Min radius (nm)", min_value=0.01, value=1.0),
            sb.number_input("Max radius (nm)", min_value=0.01, value=20.0),
        ),
        "threshold": sb.slider("Threshold", 0.005, 0.5, 0.05, step=0.005),
    }
    if detector == "log":
        detect_kw["split_touching"] = sb.checkbox("Split touching particles")
    mask_kw = _mask_controls(sb)
    sb.subheader("Output")
    thickness = sb.number_input("Section thickness (nm; 0 = no 3D estimate)", 0.0)
    bin_nm = sb.number_input("Depth bin (nm)", min_value=0.1, value=20.0)
    notes_dir = sb.text_input("Notes folder", "ferrospot-notes")

    mask, p, nb, stats, profile = _analyse(
        image, nm_per_px, mask_kw, detect_kw, thickness or None, bin_nm
    )
    stats = {k: v for k, v in asdict(stats).items() if v is not None}
    params = {"nm_per_px": nm_per_px, **detect_kw, **mask_kw}
    if thickness:
        params["section_thickness_nm"] = thickness

    cols = st.columns(6)
    cols[0].metric("Particles", stats["count"])
    for col, label, key in zip(
        cols[1:4],
        ["Diameter (nm)", "NN centre (nm)", "Gap (nm)"],
        ["d", "nn_centre", "nn_gap"],
    ):
        if f"{key}_mean_nm" in stats:
            col.metric(
                label, f"{stats[f'{key}_mean_nm']:.2f} ± {stats[f'{key}_std_nm']:.2f}"
            )
    if "overlap_fraction" in stats:
        cols[4].metric("Overlapping", f"{100 * stats['overlap_fraction']:.0f} %")
    if "clark_evans" in stats:
        cols[5].metric("Clark-Evans", f"{stats['clark_evans']:.2f}")

    left, right = st.columns([3, 1])
    left.plotly_chart(_figure(image, p, mask, nb))
    left.caption(
        "Circles: half-maximum size. Markers: red counted in the spacing "
        "statistics, orange dropped by the border correction. Cyan: mask edge."
    )
    right.dataframe(
        {"quantity": list(stats), "value": list(stats.values())}, hide_index=True
    )
    right.download_button(
        "Particles CSV",
        p.to_csv(**{k: nb[k] for k in ("nn_nm", "gap_nm", "edge_nm", "counted")}),
        file_name=f"{Path(name).stem}_particles.csv",
        mime="text/csv",
    )
    cfg = {
        "name": Path(name).stem,
        "source": {"path": name},
        "scale": {"nm_per_px": nm_per_px},
        "mask": mask_kw,
        "detect": detect_kw,
    }
    right.download_button(
        "Config JSON",
        json.dumps(cfg, indent=2),
        file_name=f"{Path(name).stem}.json",
        help="Starting point for `ferrospot run`; add the notes fields.",
    )

    if profile:
        st.subheader("Depth profile")
        fig = px.line(
            profile,
            x="depth_lo_nm",
            y=["nn_centre_mean_nm", "nn_gap_mean_nm", "d_mean_nm"],
            markers=True,
            labels={"depth_lo_nm": "depth below vacuum (bin start, nm)"},
        )
        st.plotly_chart(fig)
        st.dataframe(profile, hide_index=True)

    st.subheader(f"Notes on {name}")
    comment = st.text_area("Comment (saved with the current parameters and stats)")
    if st.button("Save note"):
        path = ferrospot.save_note(notes_dir, name, comment, params, stats)
        st.success(f"Saved to {path}")
    for f in sorted(Path(notes_dir).glob(f"{Path(name).stem}_*.md"), reverse=True):
        with st.expander(f.name):
            st.markdown(f.read_text())


if __name__ == "__main__":
    app()
