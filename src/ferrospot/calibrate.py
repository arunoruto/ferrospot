"""Forward-model calibration of 2D spacings to true 3D spacings.

A TEM image projects a section 50 to 100 nm thick, so particles at
different depths overlap and the measured 2D nearest-neighbour distance
understates the 3D one; detection pulls the other way (merged and missed
particles). Rather than correct either analytically, :func:`calibrate_panel`
simulates sections with known 3D packings in the image conditions of one
panel (pixel size, mode, blur, noise matched to the real image, the same
detector, radius range and threshold), runs the same pipeline on them and
records measured 2D against true 3D statistics; :func:`estimate_3d` then
reads the true 3D state off at the panel's measured 2D NN distance.

Per panel:

1. Sizes: lognormal, the geometric spread from the measured diameters,
   the median scaled until the simulated half-maximum diameter matches
   the measured one, cut off below the detector's smallest radius: the
   calibrated population is the particles the panel's settings count,
   not any smaller ones.
2. Noise: the pixel-scale residual (:func:`ferrospot.simulate.noise_level`)
   of the simulated image is matched to that of the real one.
3. States, at each section thickness and seed: saturated packings with a
   minimum gap of ``gaps_r`` radii ("regular") and random thinnings of a
   touching-allowed saturated packing to the same densities and below
   ("random").
4. Inverse (:func:`estimate_3d`): per thickness and family, the seed
   means form a curve of measured NN distance against true 3D gap; the
   states within 10 % of the panel's value (its sensitivity to the
   threshold) plus the seed scatter are consistent with it. The range
   over thicknesses and both families is the uncertainty.

Options, off by default: ``alphas`` adds a family of saturated packings
with radius-scaled exclusion (``alpha (r_i + r_j)``, :func:`ferrospot.
simulate.section`), and ``alpha`` to the estimate; ``bragg`` (bright
field only) simulates every state also with random Bragg visibility
(families ``<family>+bragg``), a nuisance folded into the range; ``depth``
a number-density profile along z. ``estimate_3d(..., match="nndd")``
matches the whole 2D NN-distance distribution instead of its mean.
"""

import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from ferrospot import simulate as sim
from ferrospot.config import analyse
from ferrospot.detect import detect, normalise
from ferrospot.stats import _pixels, neighbours, spacing_stats

# Measured statistics recorded per simulated image.
MEASURED = [
    "count", "n_a_per_um2", "d_mean_nm", "nn_centre_mean_nm", "nn_gap_mean_nm",
    "overlap_fraction", "clark_evans",
]  # fmt: skip
# True 3D quantities estimated from the measured NN distance.
TRUE = ["gap3d_nm", "gap3d_r", "nn3d_nm", "within_1r", "fill", "n_v_per_um3"]
# Reported by estimate_3d: true 3D state, detection completeness, and the
# simulated 2D statistics there (to compare with the measured ones).
ESTIMATED = [
    *TRUE, "completeness", "count", "n_a_per_um2", "d_mean_nm",
    "overlap_fraction", "clark_evans", "contrast_median", "alpha",
]  # fmt: skip
# Relative uncertainty of a measured mean NN distance: its change for a
# threshold change of about 0.005 (measured on real panels).
REL_ERR = 0.1
# Quantile levels recording a 2D NN-distance distribution.
QUANTILES = np.linspace(0.05, 0.95, 19)


def nn_sample(p, mask) -> np.ndarray:
    """Border-corrected 2D NN distances of the particles in the mask (the
    values whose mean is ``nn_centre_mean_nm``)."""
    p = p.subset(mask[_pixels(p, mask.shape, p.nm_per_px)])
    nb = neighbours(p, mask)
    return nb["nn_nm"][nb["counted"]]


def nn_quantiles(sample) -> list | None:
    """The :data:`QUANTILES` of a NN-distance sample (None below 5 values)."""
    sample = np.asarray(sample, float)
    return None if len(sample) < 5 else np.quantile(sample, QUANTILES).tolist()


def nndd(q, q_ref) -> float:
    """Distance of two NN-distance distributions given as quantiles.

    The Wasserstein-1 distance over the central 90 % (mean absolute
    difference of the quantiles) relative to the reference's mean: 0 for
    identical distributions, 0.1 for one scaled by 10 %.
    """
    q, q_ref = np.asarray(q, float), np.asarray(q_ref, float)
    return float(np.mean(np.abs(q - q_ref)) / np.mean(q_ref))


@dataclass
class Calibration:
    """Simulated images of one panel: settings and one row per image."""

    name: str
    nm_per_px: float
    mode: str
    detect_kw: dict
    psf_fwhm_nm: float
    noise: float
    gmean_nm: float
    gstd: float
    thickness_sd: float
    size_nm: float
    rmin_nm: float
    thicknesses_nm: list
    measured: dict = field(default_factory=dict)
    pilot: list = field(default_factory=list)
    rows: list = field(default_factory=list)
    # alphas, bragg, depth (see calibrate_panel); empty: the default grid.
    options: dict = field(default_factory=dict)

    def save(self, path):
        """Write as JSON."""
        Path(path).write_text(json.dumps(asdict(self), indent=1) + "\n")

    @classmethod
    def load(cls, path):
        """Read a calibration written by :meth:`save`."""
        return cls(**json.loads(Path(path).read_text()))


def run_one(cal: Calibration, sec: sim.Section, seed: int, bragg=None) -> dict:
    """Simulate, detect and measure one section; true and measured stats."""
    img = sim.simulate_image(
        sec, cal.nm_per_px, mode=cal.mode, psf_fwhm_nm=cal.psf_fwhm_nm,
        noise=cal.noise, thickness_sd=cal.thickness_sd, bragg=bragg, seed=seed,
    )  # fmt: skip
    mask = np.ones(img.shape, bool)
    p = detect(img, cal.nm_per_px, mask=mask, **cal.detect_kw)
    st = spacing_stats(p, mask=mask)
    row = {k: getattr(st, k) for k in MEASURED} | sim.truth(sec)
    row["nn_q"] = nn_quantiles(nn_sample(p, mask)) if len(p) >= 2 else None
    row["contrast_median"] = float(np.median(p.contrast)) if len(p) else None
    row["completeness"], row["false_fraction"] = _match(
        sec, p, img.shape[0] * cal.nm_per_px
    )
    return row


def _match(sec, p, size):
    """Share of the spheres in the section that a detection sits on, and
    share of detections on no sphere (within its radius or one pixel)."""
    vis = sec.visible()
    if not len(p) or not vis.any():
        return 0.0, None
    xy = np.column_stack([sec.x[vis], sec.y[vis]]) % size
    dist, idx = cKDTree(xy, boxsize=size).query(
        np.column_stack([p.x_nm, p.y_nm]) % size
    )
    hit = dist <= np.maximum(sec.r[vis][idx], p.nm_per_px)
    return len(np.unique(idx[hit])) / vis.sum(), float(1 - hit.mean())


def _states(args):
    """All states for one thickness and seed (a worker task)."""
    cal, t, seed, gaps_r, sparse = args
    opt = cal.options
    radii = sim.radius_spec(cal.gmean_nm, cal.gstd, cal.rmin_nm)
    r = cal.gmean_nm

    def section(**kw):
        return sim.section(
            cal.size_nm, t, radii, depth=opt.get("depth"), seed=seed, **kw
        )

    dense = section()
    secs = []
    n_last = None
    for g in gaps_r:
        try:
            reg = section(gap_nm=g * r)
        except ValueError:  # gap too large for the periodic box
            continue
        n_last = len(reg.r)
        secs.append(("regular", g, reg))
        keep = len(reg.r) / len(dense.r)
        rnd = dense.subset(np.random.default_rng(seed).random(len(dense.r)) < keep)
        secs.append(("random", keep, rnd))
    for f in sparse:
        keep = f * (n_last or len(dense.r)) / len(dense.r)
        rnd = dense.subset(np.random.default_rng(seed + 1).random(len(dense.r)) < keep)
        secs.append(("random", keep, rnd))
    for a in opt.get("alphas", ()):
        try:
            secs.append(("alpha", a, section(packing="alpha", alpha=a)))
        except ValueError:
            continue
    braggs = [None] + ([opt["bragg"]] if opt.get("bragg") and cal.mode == "bf" else [])
    rows = []
    for fam, param, sec in secs:
        for b in braggs:
            row = {"family": fam + ("+bragg" if b else ""), "param": param}
            row |= run_one(cal, sec, seed, bragg=b)
            if fam == "alpha":
                row["alpha"] = param
            rows.append(row)
    for row in rows:
        row.update(thickness_nm=t, seed=seed)
    return rows


def match_noise(cal: Calibration, target: float, sec: sim.Section) -> float:
    """White-noise level that gives a simulated image the target residual.

    Noise adds in quadrature to the structure of the noise-free image;
    the residual filter's gain for white noise is measured on pure noise.
    """
    clean = sim.simulate_image(
        sec, cal.nm_per_px, mode=cal.mode, psf_fwhm_nm=cal.psf_fwhm_nm, noise=0.0,
        thickness_sd=cal.thickness_sd,
    )  # fmt: skip
    e0 = sim.noise_level(clean)
    gain = sim.noise_level(np.random.default_rng(0).normal(size=(256, 256)))
    return float(np.sqrt(max(target**2 - e0**2, 0.0)) / gain)


def calibrate_panel(
    cfg: dict,
    result: dict | None = None,
    *,
    base=".",
    thicknesses_nm=None,
    gaps_r=(0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16),
    sparse=(0.5, 0.25),
    seeds=(0, 1, 2),
    image_px: int = 300,
    psf_fwhm_nm: float | None = None,
    workers: int = 1,
    alphas=(),
    bragg=None,
    depth: dict | None = None,
) -> Calibration:
    """Simulate a panel's imaging conditions over a grid of 3D states.

    Parameters
    ----------
    cfg : dict
        Panel config (see :mod:`ferrospot.config`); uses ``detect`` and
        ``section_thickness_nm``.
    result : dict, optional
        Output of :func:`ferrospot.analyse` for it (needs ``image``,
        ``mask``, ``nm_per_px`` and ``stats``); run if not given.
    base : str or Path
        Folder the config's relative paths are relative to.
    thicknesses_nm : sequence of float, optional
        Section thicknesses to simulate. Default: 0.8, 1 and 1.2 times
        ``section_thickness_nm``, or 50, 75 and 100 nm without it.
    gaps_r : sequence of float
        Minimum gaps of the saturated ("regular") states, in median radii.
    sparse : sequence of float
        Extra random states at these fractions of the sparsest regular
        density.
    seeds : sequence of int
        One packing and noise realisation per seed and state.
    image_px : int
        Simulated image edge in pixels, enlarged up to twice that when
        the largest gap does not fit the periodic box; gaps that still
        do not fit are left out.
    psf_fwhm_nm : float, optional
        Blur FWHM. Default: 2 pixels.
    workers : int
        Processes.
    alphas : sequence of float
        Extra family "alpha": saturated packings with exclusion
        ``alpha (r_i + r_j)``. Default: none.
    bragg : dict or True, optional
        Bright-field panels only: every state also with random Bragg
        visibility (:func:`ferrospot.simulate.bragg_factors`), as families
        ``<family>+bragg``. The pilot runs without it.
    depth : dict, optional
        Number-density profile along z (:func:`ferrospot.simulate.section`)
        for every state and the pilot.

    Returns
    -------
    Calibration
    """
    if result is None:
        result = analyse(cfg, base)
    s, st = result["nm_per_px"], result["stats"]
    kw = dict(cfg.get("detect", {}))
    bright_field = kw.get("bright_field", False)
    t0 = cfg.get("section_thickness_nm")
    if thicknesses_nm is None:
        thicknesses_nm = [0.8 * t0, t0, 1.2 * t0] if t0 else [50.0, 75.0, 100.0]
    d, sd = st.d_mean_nm, st.d_std_nm or 0.0
    gstd = float(np.exp(np.sqrt(np.log1p((sd / d) ** 2))))
    # Half-maximum outline of a sharp sphere: sqrt(3)/2 of its diameter.
    gmean = d / np.sqrt(3) / np.exp(0.5 * np.log(gstd) ** 2)
    rmax = gmean * gstd**2
    size = max(image_px * s, 2.2 * (2 * rmax + max(gaps_r) * gmean))
    # ponytail: larger gaps than fit a 2 * image_px box are skipped, not simulated.
    size = min(size, 2 * image_px * s)
    size = np.ceil(size / s) * s
    cal = Calibration(
        name=cfg.get("name", ""), nm_per_px=s, mode="bf" if bright_field else "haadf",
        detect_kw=kw, psf_fwhm_nm=psf_fwhm_nm or 2 * s, noise=0.0, gmean_nm=gmean,
        gstd=gstd, thickness_sd=0.05, size_nm=float(size), thicknesses_nm=list(thicknesses_nm),
        rmin_nm=float(kw.get("radius_nm", (0.0,))[0]),
        measured={k: getattr(st, k) for k in MEASURED},
    )  # fmt: skip
    if bragg is True:
        bragg = dict(sim.BRAGG)
    opts = {"alphas": list(alphas), "bragg": bragg, "depth": depth}
    cal.options = {k: v for k, v in opts.items() if v}
    real = normalise(result["image"], result["mask"], bright_field)
    target = sim.noise_level(real, result["mask"])
    cal.measured["noise_residual"] = target
    if "particles" in result:
        cal.measured["contrast_median"] = float(np.median(result["particles"].contrast))
        cal.measured["nn_q"] = nn_quantiles(
            nn_sample(result["particles"], result["mask"])
        )

    # Pilot: a random state, starting at the density a fully detected
    # Poisson pattern would have (N_A / (t + d)). Match, in turn, the pixel
    # noise, the median particle contrast (with the thickness variation,
    # which sets how much of the stretched grey range the particles get),
    # the mean diameter (with the median radius) and the mean NN distance
    # (with the density), so the first three are matched near the state
    # the panel will be assigned.
    t_mid = thicknesses_nm[len(thicknesses_nm) // 2]
    n_v = st.n_a_per_um2 * 1e-6 / (t_mid + d)
    c_real, nn_real = cal.measured.get("contrast_median"), st.nn_centre_mean_nm
    for _ in range(8):
        radii = sim.radius_spec(cal.gmean_nm, cal.gstd, cal.rmin_nm)
        dense = sim.section(cal.size_nm, t_mid, radii, depth=depth, seed=99)
        n_dense = sim.truth(dense)["n_v_per_um3"] * 1e-9
        n_v = min(n_v, n_dense)
        pilot = dense.subset(
            np.random.default_rng(99).random(len(dense.r)) < n_v / n_dense
        )
        cal.noise = match_noise(cal, target, pilot)
        row = run_one(cal, pilot, 99)
        if not row["d_mean_nm"] or not row["nn_centre_mean_nm"]:
            n_v *= 3
            continue
        cal.pilot.append(
            {"gmean_nm": cal.gmean_nm, "noise": cal.noise,
             "thickness_sd": cal.thickness_sd, "n_v_per_um3": n_v * 1e9}
            | {k: row[k] for k in ("d_mean_nm", "nn_centre_mean_nm", "contrast_median")}
        )  # fmt: skip
        # Bounded: in very noisy panels noise blobs set the measured size.
        f = np.clip(d / row["d_mean_nm"], 0.7, 1.4)
        cal.gmean_nm = float(np.clip(cal.gmean_nm * f, 0.6 * gmean, 1.6 * gmean))
        if c_real:
            # c = A / (A + range): scale the range by the change in 1/c - 1.
            f = (1 / c_real - 1) / (1 / row["contrast_median"] - 1)
            cal.thickness_sd = float(np.clip(cal.thickness_sd * f, 0.0, 0.3))
        n_v *= float(np.clip((row["nn_centre_mean_nm"] / nn_real) ** 3, 0.3, 3.0))

    tasks = [(cal, t, seed, gaps_r, sparse) for t in thicknesses_nm for seed in seeds]
    if workers > 1:
        with ProcessPoolExecutor(workers) as ex:
            parts = list(ex.map(_states, tasks))
    else:
        parts = [_states(a) for a in tasks]
    cal.rows = [row for part in parts for row in part]
    return cal


def curves(cal: Calibration, key="nn_centre_mean_nm") -> list[dict]:
    """Seed means per state: one curve per thickness and family.

    Returns
    -------
    list of dict
        ``thickness_nm``, ``family`` and, per state sorted by the true 3D
        gap, arrays of the seed means of every recorded number and
        ``ln_sd``, the seed standard deviation of ``ln(key)``.
    """
    out = []
    for t in cal.thicknesses_nm:
        for fam in sorted({r["family"] for r in cal.rows}):
            rows = [
                r for r in cal.rows if r["thickness_nm"] == t and r["family"] == fam
            ]
            states = sorted({round(r["param"], 9) for r in rows})
            pts = []
            for st in states:
                grp = [r for r in rows if round(r["param"], 9) == st and r.get(key)]
                if not grp:
                    continue
                num = {
                    k: np.mean([np.nan if r.get(k) is None else r[k] for r in grp])
                    for k in grp[0]
                    if k not in ("family", "seed", "thickness_nm", "nn_q")
                }
                if "nn_q" in grp[0]:
                    qs = [r["nn_q"] for r in grp if r.get("nn_q")]
                    num["nn_q"] = (
                        np.mean(qs, axis=0) if qs else np.full(len(QUANTILES), np.nan)
                    )
                num["ln_sd"] = float(np.std(np.log([r[key] for r in grp])))
                if np.isfinite(num["gap3d_r"]):
                    pts.append(num)
            pts.sort(key=lambda p: p["gap3d_r"])
            if len(pts) >= 2:
                c = {k: np.array([p[k] for p in pts]) for k in pts[0]}
                out.append({"thickness_nm": t, "family": fam, **c})
    return out


def estimate_3d(
    measured, cal: Calibration, key="nn_centre_mean_nm", rel_err=REL_ERR,
    match="mean",
) -> dict:  # fmt: skip
    """True 3D state consistent with a measured 2D statistic.

    Each curve of :func:`curves` is interpolated in ``ln`` of the true 3D
    gap; its states whose simulated ``ln(key)`` lies within ``rel_err``
    (added in quadrature to the seed scatter) of the measured value are
    consistent. The estimate on a curve is where it comes closest to the
    measured value; ``_lo`` / ``_hi`` span all consistent states of all
    thicknesses and both families. Where the curve is flat (detection
    saturated) the range opens up to the end of the simulated grid.

    ``match="nndd"`` compares the whole NN-distance distribution instead
    (:func:`nndd` of the simulated and the measured quantiles,
    ``cal.measured["nn_q"]``): the estimate is the state of least
    distance, consistent are those within ``rel_err`` (plus the seed
    scatter of ``ln(key)``) of that least distance, reported as ``nndd``.

    Parameters
    ----------
    measured : SpacingStats or dict or float
        The panel's statistics, or the value of ``key`` itself.
    cal : Calibration
        From :func:`calibrate_panel` for the same panel.
    key : str
        Measured statistic to invert (default the mean NN distance).
    rel_err : float
        Relative uncertainty of the measured value.
    match : {"mean", "nndd"}
        Match ``key``, or the NN-distance distribution.

    Returns
    -------
    dict
        ``<q>``, ``<q>_lo``, ``<q>_hi`` for every ``q`` in ``ESTIMATED``:
        the mean of the two families' estimates at the middle thickness,
        and the consistent range; ``in_range`` (the measured value lies
        within every curve's simulated range) and ``by_curve`` (estimate
        and range per thickness and family); with ``match="nndd"`` also
        ``nndd``, the least distance at the middle thickness.
    """
    if match not in ("mean", "nndd"):
        raise ValueError(f"match must be 'mean' or 'nndd', not {match!r}")
    q_ref = cal.measured.get("nn_q") if match == "nndd" else None
    if match == "nndd" and not q_ref:
        raise ValueError("no measured NN quantiles: recalibrate (cal.measured['nn_q'])")
    if isinstance(measured, (int, float)):
        x0 = float(measured)
    elif isinstance(measured, dict):
        x0 = measured[key]
    else:
        x0 = getattr(measured, key)
    lx = np.log(x0)
    by_curve = []
    for c in curves(cal, key):
        lg = np.log(c["gap3d_r"])
        fine = np.linspace(lg[0], lg[-1], 400)
        y = {
            q: np.interp(fine, lg, c[q])
            for q in ESTIMATED
            if q in c and np.isfinite(c[q]).all()
        }
        m = np.interp(fine, lg, np.log(c[key]))
        tol = np.hypot(rel_err, np.interp(fine, lg, c["ln_sd"]))
        est = {"thickness_nm": c["thickness_nm"], "family": c["family"]}
        if q_ref:
            dist = np.array([nndd(q, q_ref) for q in c.get("nn_q", [])])
            fin = np.isfinite(dist)
            if fin.sum() < 2:
                continue
            d = np.interp(fine, lg[fin], dist[fin])
            est["nndd"] = float(d.min())
            d = d - d.min()
        else:
            d = np.abs(m - lx)
        best = np.argmin(d)
        ok = d <= tol
        ok[best] = True
        est["in_range"] = bool(m.min() <= lx <= m.max())
        for q, v in y.items():
            est[q], est[q + "_lo"], est[q + "_hi"] = (
                float(v[best]), float(v[ok].min()), float(v[ok].max())
            )  # fmt: skip
        by_curve.append(est)
    out = {"measured": x0, "key": key, "match": match, "by_curve": by_curve}
    t_mid = cal.thicknesses_nm[len(cal.thicknesses_nm) // 2]
    for q in [*ESTIMATED, "nndd"]:
        mid = [e[q] for e in by_curve if e["thickness_nm"] == t_mid and q in e]
        if mid:
            out[q] = float(np.mean(mid))
            if q == "nndd":
                continue
            out[q + "_lo"] = min(e[q + "_lo"] for e in by_curve if q in e)
            out[q + "_hi"] = max(e[q + "_hi"] for e in by_curve if q in e)
    out["in_range"] = all(e["in_range"] for e in by_curve)
    return out


def simulate_state(cal: Calibration, n_v_per_um3: float, thickness_nm: float, seed=0):
    """Random section at a given number density, imaged like the panel.

    Returns
    -------
    sec : ferrospot.simulate.Section
    image : numpy.ndarray
    """
    radii = sim.radius_spec(cal.gmean_nm, cal.gstd, cal.rmin_nm)
    depth = cal.options.get("depth")
    dense = sim.section(cal.size_nm, thickness_nm, radii, depth=depth, seed=seed)
    n_dense = sim.truth(dense)["n_v_per_um3"]
    keep = np.random.default_rng(seed).random(len(dense.r)) < n_v_per_um3 / n_dense
    sec = dense.subset(keep)
    img = sim.simulate_image(
        sec, cal.nm_per_px, mode=cal.mode, psf_fwhm_nm=cal.psf_fwhm_nm,
        noise=cal.noise, thickness_sd=cal.thickness_sd, seed=seed,
    )  # fmt: skip
    return sec, img


def table(entries) -> str:
    """Markdown table of ``(name, calibration, estimate)`` triples.

    The estimate inverts the mean NN distance; the column "from N_A"
    inverts the measured count per area instead. Where the two disagree,
    the panel is less uniform (or more clustered) than either family of
    simulated states.
    """

    def rng(e, q, f="{:.1f}", scale=1.0):
        v = [e.get(k) for k in (q, q + "_lo", q + "_hi")]
        if v[0] is None:
            return "–"
        a, lo, hi = (scale * x for x in v)
        return f"{f.format(a)} ({f.format(lo)}–{f.format(hi)})"

    head = (
        "| panel | t (nm) | 2D NN (nm) | 3D gap (nm) | 3D gap / r | gap < r |"
        " 3D NN (nm) | fill (%) | n_V (µm⁻³) | completeness | sim / real: N_A, overlap, CE |"
        " 3D gap / r from N_A |"
    )
    with_q = any(cal.measured.get("nn_q") for _, cal, _ in entries)
    with_a = any(e.get("alpha") is not None for *_, e in entries)
    head += " 3D gap / r from NNDD (misfit) |" * with_q + " alpha |" * with_a
    lines = [head, "|" + "---|" * (12 + with_q + with_a)]
    for name, cal, e in entries:
        t = cal.thicknesses_nm
        m = cal.measured
        check = ", ".join(
            "–" if e.get(k) is None or m.get(k) is None else f"{e[k] / m[k]:.2f}"
            for k in ("n_a_per_um2", "overlap_fraction", "clark_evans")
        )
        flag = "" if e["in_range"] else " ⚠"
        lines.append(
            f"| {name}{flag} | {t[0]:.0f}–{t[-1]:.0f} | {e['measured']:.2f} "
            f"| {rng(e, 'gap3d_nm')} | {rng(e, 'gap3d_r')} "
            f"| {rng(e, 'within_1r', '{:.0f}', 100)} % | {rng(e, 'nn3d_nm')} "
            f"| {rng(e, 'fill', '{:.2f}', 100)} | {rng(e, 'n_v_per_um3', '{:.2g}')} "
            f"| {rng(e, 'completeness', '{:.0f}', 100)} % | {check} "
            f"| {rng(estimate_3d(m, cal, 'n_a_per_um2'), 'gap3d_r')} |"
        )
        if with_q:
            q = estimate_3d(m, cal, match="nndd") if m.get("nn_q") else {}
            mis = f" ({q['nndd']:.2f})" if "nndd" in q else ""
            lines[-1] += f" {rng(q, 'gap3d_r')}{mis} |"
        if with_a:
            lines[-1] += f" {rng(e, 'alpha', '{:.2f}')} |"
    return "\n".join(lines)


def plot_calibration(cal: Calibration, est: dict | None = None, axes=None):
    """Calibration curves: measured 2D NN distance, overlap and completeness
    against the true 3D gap, per thickness (colour) and family (line
    style), with the panel's value and the inferred range.

    Needs matplotlib (``simulate`` extra). Returns the figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    if axes is None:
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    fig = axes[0].figure
    colours = plt.cm.viridis(np.linspace(0, 0.85, len(cal.thicknesses_nm)))
    panels = [
        ("nn_centre_mean_nm", "measured 2D NN distance (nm)"),
        ("overlap_fraction", "measured overlap in projection"),
        ("completeness", "detection completeness"),
    ]
    for c in curves(cal):
        col = colours[cal.thicknesses_nm.index(c["thickness_nm"])]
        ls = {"random": "-", "regular": "--"}.get(c["family"].split("+")[0], ":")
        for ax, (k, _) in zip(axes, panels):
            label = f"t = {c['thickness_nm']:.0f} nm, {c['family']}"
            ax.plot(c["gap3d_r"], c[k], ls, marker="o", ms=3, color=col, label=label)
    m = cal.measured
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for ax, (k, ylabel) in zip(axes, panels):
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
        ax.xaxis.set_major_formatter(plain)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("true 3D NN surface gap / radius")
        ax.set_ylabel(ylabel)
        if m.get(k) is not None:
            ax.axhline(m[k], color="crimson", lw=1)
        if est and est.get("gap3d_r") is not None:
            ax.axvspan(
                est["gap3d_r_lo"], est["gap3d_r_hi"], color="crimson", alpha=0.12
            )
            ax.axvline(est["gap3d_r"], color="crimson", lw=1, ls=":")
    x0 = m.get("nn_centre_mean_nm")
    if x0:
        axes[0].axhspan(
            x0 * (1 - REL_ERR), x0 * (1 + REL_ERR), color="crimson", alpha=0.12
        )
    axes[0].set_yscale("log")
    axes[0].yaxis.set_major_locator(LogLocator(subs=(1, 2, 3, 5)))
    axes[0].yaxis.set_major_formatter(plain)
    axes[0].yaxis.set_minor_formatter(NullFormatter())
    axes[0].legend(fontsize=7)
    fig.suptitle(
        f"{cal.name}: {cal.mode.upper()}, {cal.nm_per_px:.3f} nm/px, "
        f"PSF {cal.psf_fwhm_nm:.2f} nm; red: measured 2D value and inferred 3D range",
        fontsize=9,
    )
    return fig


def plot_strip(cal: Calibration, gaps_r=(0.5, 1, 2, 4, 8), thickness_nm=None, seed=0):
    """Simulated sections at increasing 3D gap (saturated packings), with
    the detections and their measured 2D NN distance. Returns the figure."""
    import matplotlib.pyplot as plt

    from ferrospot.image import overlay

    t = thickness_nm or cal.thicknesses_nm[len(cal.thicknesses_nm) // 2]
    radii = sim.radius_spec(cal.gmean_nm, cal.gstd, cal.rmin_nm)
    fig, axes = plt.subplots(
        1, len(gaps_r), figsize=(2.9 * len(gaps_r), 3.4), layout="constrained"
    )
    for ax, g in zip(axes, gaps_r):
        sec = sim.section(cal.size_nm, t, radii, gap_nm=g * cal.gmean_nm, seed=seed)
        img = sim.simulate_image(
            sec, cal.nm_per_px, mode=cal.mode, psf_fwhm_nm=cal.psf_fwhm_nm,
            noise=cal.noise, thickness_sd=cal.thickness_sd, seed=seed,
        )  # fmt: skip
        mask = np.ones(img.shape, bool)
        p = detect(img, cal.nm_per_px, mask=mask, **cal.detect_kw)
        st, tr = spacing_stats(p, mask=mask), sim.truth(sec)
        ax.imshow(overlay(img, p))
        ax.set_axis_off()
        ax.set_title(
            f"3D gap {tr['gap3d_r']:.1f} r ({tr['gap3d_nm']:.1f} nm)\n"
            f"2D NN {st.nn_centre_mean_nm:.1f} nm, overlap {100 * st.overlap_fraction:.0f} %",
            fontsize=8,
        )
    fig.suptitle(
        f"Simulated {t:.0f} nm sections, {cal.size_nm:.0f} nm wide, imaged like {cal.name}",
        fontsize=9,
    )
    return fig


def plot_compare(result: dict, cal: Calibration, est: dict, seed=0):
    """Real panel (a crop the size of the simulation) next to a simulated
    section at the inferred 3D density, both with detections. Returns the
    figure. For local use: published panels are copyrighted."""
    import matplotlib.pyplot as plt

    from ferrospot.image import overlay

    t = cal.thicknesses_nm[len(cal.thicknesses_nm) // 2]
    _, img = simulate_state(cal, est["n_v_per_um3"], t, seed)
    mask = np.ones(img.shape, bool)
    p = detect(img, cal.nm_per_px, mask=mask, **cal.detect_kw)
    real, rmask, rp = result["image"], result["mask"], result["particles"]
    n = min(img.shape[0], *real.shape)
    # The crop with the most mask in it.
    cover = np.cumsum(np.cumsum(np.pad(rmask, ((1, 0), (1, 0))), 0), 1)
    sums = cover[n:, n:] - cover[:-n, n:] - cover[n:, :-n] + cover[:-n, :-n]
    y0, x0 = np.unravel_index(np.argmax(sums), sums.shape)
    s = cal.nm_per_px
    inside = (
        (rp.x_nm >= x0 * s) & (rp.x_nm < (x0 + n) * s)
        & (rp.y_nm >= y0 * s) & (rp.y_nm < (y0 + n) * s)
    )  # fmt: skip
    crop_p = rp.subset(inside)
    crop_p.x_nm, crop_p.y_nm = crop_p.x_nm - x0 * s, crop_p.y_nm - y0 * s
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.3), layout="constrained")
    axes[0].imshow(
        overlay(real[y0 : y0 + n, x0 : x0 + n], crop_p, rmask[y0 : y0 + n, x0 : x0 + n])
    )
    axes[0].set_title(f"{cal.name} (crop {n * s:.0f} nm)", fontsize=9)
    sub = p.subset((p.x_nm < n * s) & (p.y_nm < n * s))
    axes[1].imshow(overlay(img[:n, :n], sub))
    axes[1].set_title(
        f"simulated: t = {t:.0f} nm, 3D gap {est['gap3d_r']:.1f} r, "
        f"fill {100 * est['fill']:.2f} %",
        fontsize=9,
    )
    for ax in axes:
        ax.set_axis_off()
    return fig
