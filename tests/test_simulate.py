"""Forward model (ferrospot.simulate) and 2D -> 3D calibration."""

import numpy as np
import pytest

import ferrospot
from ferrospot import calibrate as cb
from ferrospot import simulate as sim

polydisc = pytest.importorskip("polydisc")


def one_sphere(z, t=20.0, r=5.0, size=40.0):
    a = np.array
    return sim.Section(a([20.0]), a([20.0]), a([z]), a([r]), size, t)


def measure(img, s, **kw):
    mask = np.ones(img.shape, bool)
    p = ferrospot.detect(img, s, mask=mask, **kw)
    return ferrospot.spacing_stats(p, mask=mask)


def test_chord_profile_and_contrast_sign():
    path = sim.iron_path(one_sphere(z=10.0), 0.5)
    rho = (np.arange(path.shape[1]) + 0.5) * 0.5 - 20.0
    row = path[40]  # pixel centre at y = 20.25 nm
    expected = 2 * np.sqrt(np.clip(25.0 - rho**2 - 0.25**2, 0, None))
    np.testing.assert_allclose(row, expected, atol=1e-9)

    # Cut by the section surface: chord from z = 0 to z + sqrt(r^2 - rho^2).
    cut = sim.iron_path(one_sphere(z=2.0), 0.5)
    assert cut[40, 40] == pytest.approx(2.0 + np.sqrt(25.0 - 2 * 0.25**2))

    for mode, sign in [("haadf", 1), ("bf", -1)]:
        img = sim.simulate_image(one_sphere(z=10.0), 0.5, mode=mode, noise=0.0)
        assert sign * (img[40, 40] - img[5, 5]) > 0.5


def test_sparse_monolayer_measures_true_spacing():
    # One layer of well-separated spheres: no overlap, 2D NN = 3D NN.
    r, size = 2.0, 120.0
    pack = polydisc.sample(
        polydisc.Box([size, size], periodic=(True, True)), r, gap=6 * r, seed=3
    )
    (x, y), radii = pack.points.T, pack.radii
    sec = sim.Section(x, y, np.full(len(x), r), radii, size, 2 * r)
    img = sim.simulate_image(sec, 0.5, noise=0.01)

    st = measure(img, 0.5, radius_nm=(1.0, 4.0), threshold=0.05)

    assert st.count == len(x)
    assert st.overlap_fraction == 0
    assert st.nn_centre_mean_nm == pytest.approx(sim.truth(sec)["nn3d_nm"], rel=0.05)


def test_overlap_grows_with_thickness():
    kw = {"radius_nm": (1.0, 4.0), "threshold": 0.05}
    overlap = []
    for t in (8.0, 80.0):
        sec = sim.section(60.0, t, 1.5, gap_nm=6.0, seed=1)
        st = measure(sim.simulate_image(sec, 0.5, noise=0.01), 0.5, **kw)
        overlap.append(st.overlap_fraction)
    assert overlap[1] > overlap[0] + 0.1


def test_deterministic():
    a = sim.section(40.0, 30.0, sim.radius_spec(1.5, 1.2), gap_nm=1.0, seed=5)
    b = sim.section(40.0, 30.0, sim.radius_spec(1.5, 1.2), gap_nm=1.0, seed=5)
    np.testing.assert_array_equal(a.x, b.x)
    img = [sim.simulate_image(a, 0.5, thickness_sd=0.05, seed=s) for s in (0, 0, 1)]
    np.testing.assert_array_equal(img[0], img[1])
    assert not np.array_equal(img[0], img[2])


def test_calibration_recovers_held_out_state():
    s, t = 0.5, 60.0
    cfg = {
        "name": "synthetic",
        "detect": {"detector": "log", "radius_nm": (1.0, 4.0), "threshold": 0.05},
        "section_thickness_nm": t,
    }
    held = sim.section(60.0, t, sim.radius_spec(2.0, 1.15), gap_nm=3 * 2.0, seed=42)
    img = sim.simulate_image(held, s, noise=0.03, seed=42)
    mask = np.ones(img.shape, bool)
    p = ferrospot.detect(img, s, mask=mask, **cfg["detect"])
    result = {
        "image": img, "mask": mask, "nm_per_px": s, "particles": p,
        "stats": ferrospot.spacing_stats(p, mask=mask),
    }  # fmt: skip

    cal = cb.calibrate_panel(
        cfg, result, thicknesses_nm=[t], gaps_r=(0.5, 1, 2, 4, 8),
        sparse=(0.5,), seeds=(0, 1), image_px=120,
    )  # fmt: skip
    est = cb.estimate_3d(result["stats"], cal)
    by_nndd = cb.estimate_3d(result["stats"], cal, match="nndd")

    true = sim.truth(held)
    assert by_nndd["gap3d_r_lo"] <= true["gap3d_r"] <= by_nndd["gap3d_r_hi"]
    assert est["in_range"]
    assert est["gap3d_r_lo"] <= true["gap3d_r"] <= est["gap3d_r_hi"]
    assert est["fill_lo"] <= true["fill"] <= est["fill_hi"]
    assert est["gap3d_r"] == pytest.approx(true["gap3d_r"], rel=0.5)


def monolayer(n_side=6, r=2.0, pitch=12.0):
    g = (np.arange(n_side) + 0.5) * pitch
    x, y = (v.ravel() for v in np.meshgrid(g, g))
    return sim.Section(
        x, y, np.full(x.size, r), np.full(x.size, r), n_side * pitch, 2 * r
    )


def test_bragg_spreads_bright_field_contrast_and_is_seeded():
    sec, s = monolayer(), 0.5
    ix = (sec.x / s).astype(int)

    def centres(**kw):
        img = sim.simulate_image(sec, s, noise=0.0, **kw)
        return img[ix, (sec.y / s).astype(int)]

    plain = centres(mode="bf")
    a, b = (
        centres(mode="bf", bragg=True, seed=3),
        centres(mode="bf", bragg=True, seed=3),
    )
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, centres(mode="bf", bragg=True, seed=4))
    assert np.std(a) > 5 * np.std(plain) + 0.02  # strong / mid / weak particles
    # HAADF: no effect by default, a small one on request.
    h = centres(mode="haadf")
    np.testing.assert_array_equal(h, centres(mode="haadf", bragg=True))
    assert not np.array_equal(h, centres(mode="haadf", bragg={"haadf": 0.2}))


def test_alpha_packing_respects_radius_scaled_exclusion():
    alpha, size = 1.5, 25.0
    sec = sim.section(
        size, 12.0, sim.radius_spec(1.0, 1.3), packing="alpha", alpha=alpha, seed=2
    )
    d = [np.subtract.outer(v, v) for v in (sec.x, sec.y, sec.z)]
    wrapped = (np.abs(d[0]) > size / 2) | (np.abs(d[1]) > size / 2)
    d[0] -= size * np.round(d[0] / size)
    d[1] -= size * np.round(d[1] / size)
    ratio = np.sqrt(sum(v**2 for v in d)) / np.add.outer(sec.r, sec.r)
    np.fill_diagonal(ratio, np.inf)
    assert sec.r.std() > 0
    assert ratio.min() >= alpha - 1e-9
    assert ratio[wrapped].min() < 1.1 * alpha  # saturated across x / y too
    with pytest.raises(ValueError):
        sim.section(size, 12.0, 1.0, packing="alpha", alpha=0.9)


def test_depth_profile_shapes_z_histogram():
    t = 60.0
    sec = sim.section(40.0, t, 1.0, depth={"kind": "exp", "scale_nm": 15.0}, seed=1)
    assert sec.z.min() >= 0
    n, _ = np.histogram(sec.z, bins=[0, 15, 30, 45])
    assert n[0] / n[1] == pytest.approx(np.e, rel=0.15)
    assert n[1] / n[2] == pytest.approx(np.e, rel=0.25)
    peak = sim.section(
        40.0, t, 1.0, depth={"kind": "peak", "z0_nm": 30.0, "w_nm": 6.0}, seed=1
    )
    inside = peak.z[(peak.z > 0) & (peak.z < t)]
    assert np.mean(inside) == pytest.approx(30.0, abs=1.0)
    assert np.std(inside) == pytest.approx(6.0, rel=0.15)


def test_nndd_distance():
    a = np.random.default_rng(0).gamma(4.0, 1.0, 500)
    q = cb.nn_quantiles(a)
    assert cb.nndd(q, q) == 0
    d = [cb.nndd(cb.nn_quantiles(a + f * a.mean()), q) for f in (0.05, 0.1, 0.3)]
    assert 0 < d[0] < d[1] < d[2]
    assert d[1] == pytest.approx(0.1, rel=0.05)
