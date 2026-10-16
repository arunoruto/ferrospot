import json
from datetime import UTC, datetime

import numpy as np
import pytest
from PIL import Image
from scipy import ndimage as ndi

import ferrospot
from ferrospot import cli

NM_PER_PX = 0.5
HALF_MAX = np.sqrt(3) / 2  # half-maximum radius of a projected sphere, per R


def spheres(shape, centres, radii, blur=1.0, noise=0.05, seed=0):
    """Projected spheres (intensity ~ chord length) plus blur and noise."""
    yy, xx = np.indices(shape)
    img = np.zeros(shape)
    for (x, y), r in zip(centres, radii):
        img += np.sqrt(np.clip(r**2 - (yy - y) ** 2 - (xx - x) ** 2, 0, None)) / r
    img = ndi.gaussian_filter(img, blur)
    return img + np.random.default_rng(seed).normal(0, noise, shape)


def match(p, centres):
    """Index of the detection nearest to each true centre, and its distance."""
    found = np.column_stack([p.x_nm, p.y_nm]) / p.nm_per_px
    d = np.hypot(*(found[:, None] - centres[None]).transpose(2, 0, 1))
    return d.argmin(0), d.min(0)


def grid(n, a, jitter, seed):
    """n x n jittered lattice of spacing a (px), starting at a."""
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(1, n + 1) * a, np.arange(1, n + 1) * a)
    return np.column_stack([gx.ravel(), gy.ravel()]) + rng.uniform(
        -jitter, jitter, (n * n, 2)
    )


def lattice_particles(n, a, r):
    """n x n square lattice (spacing a nm, radius r nm) at 1 nm/px."""
    gx, gy = np.meshgrid(np.arange(n) * a + a / 2, np.arange(n) * a + a / 2)
    k = n * n
    return ferrospot.Particles(
        gx.ravel(), gy.ravel(), np.full(k, r), np.ones(k), 1.0
    )  # fmt: skip


def test_log_recovers_positions_and_half_max_radii():
    centres = grid(5, 40, 3, seed=1) - 10
    radii = np.random.default_rng(1).uniform(4, 10, 25)  # px
    img = spheres((220, 220), centres, radii)

    p = ferrospot.detect(img, NM_PER_PX, radius_nm=(1.0, 6.0), threshold=0.1)

    assert len(p) == 25
    idx, dist = match(p, centres)
    assert dist.max() < 1.0
    assert np.allclose(p.r_nm[idx] / NM_PER_PX, HALF_MAX * radii, rtol=0.15)
    assert (p.contrast > 0.3).all()


@pytest.mark.parametrize("r", [4.0, 8.0, 16.0])
def test_half_max_size_of_spheres(r):
    img = spheres((100, 100), [(50, 50)], [r], blur=0.5, noise=0.02)
    p = ferrospot.detect(img, NM_PER_PX, radius_nm=(1.0, 10.0), threshold=0.05)
    assert len(p) == 1
    assert p.r_nm[0] / NM_PER_PX == pytest.approx(HALF_MAX * r, rel=0.08)


def test_dog_finds_dense_small_particles():
    # 3 px radius, centres 8 px apart: gaps of 2 px, blurred into each other.
    centres = grid(20, 8, 1, seed=2)
    img = spheres((172, 172), centres, np.full(len(centres), 3.0), noise=0.03)

    p = ferrospot.detect(
        img, NM_PER_PX, detector="dog", radius_nm=(1.0, 2.5), threshold=0.02
    )

    _, dist = match(p, centres)
    assert np.mean(dist < 1.5) > 0.95
    assert len(p) <= 1.05 * len(centres)
    assert np.median(p.r_nm / NM_PER_PX) == pytest.approx(HALF_MAX * 3, rel=0.15)


def test_bright_field_inverts():
    centres = np.array([(30.0, 30.0), (80.0, 40.0), (50.0, 80.0)])
    img = 1.0 - spheres((110, 110), centres, [6, 8, 5])  # dark particles
    kw = {"radius_nm": (1.5, 6.0), "threshold": 0.1}

    p = ferrospot.detect(img, NM_PER_PX, bright_field=True, **kw)

    assert len(p) == 3
    assert match(p, centres)[1].max() < 1.0
    # Without inversion only noise turns up, none of it on a particle.
    assert match(ferrospot.detect(img, NM_PER_PX, **kw), centres)[1].min() > 5


def test_split_touching_pair():
    r = 8
    img = spheres((60, 80), [(40 - 0.8 * r, 30), (40 + 0.8 * r, 30)], [r, r])
    kw = {"radius_nm": (2.0, 10.0), "threshold": 0.1}
    assert len(ferrospot.detect(img, NM_PER_PX, **kw)) == 1
    p = ferrospot.detect(img, NM_PER_PX, split_touching=True, **kw)
    assert len(p) == 2
    assert np.allclose(
        np.sort(p.x_nm) / NM_PER_PX, [40 - 0.8 * r, 40 + 0.8 * r], atol=1.5
    )
    # Each half-maximum region stops where the Voronoi cells meet.
    assert np.allclose(p.r_nm / NM_PER_PX, HALF_MAX * r, rtol=0.25)


def test_split_leaves_single_particles_alone():
    img = spheres((60, 60), [(30, 30)], [8])
    p = ferrospot.detect(img, NM_PER_PX, radius_nm=(2.0, 10.0), split_touching=True)
    assert len(p) == 1


def test_mask_drops_particles_outside():
    img = spheres((60, 120), [(30, 30), (90, 30)], [6, 6])
    mask = np.zeros(img.shape, bool)
    mask[:, :60] = True
    p = ferrospot.detect(img, NM_PER_PX, radius_nm=(1.0, 5.0), threshold=0.1, mask=mask)
    assert len(p) == 1 and p.x_nm[0] < 60 * NM_PER_PX


def test_mask_edge_inside_material_is_not_a_particle():
    # A depth band cuts through a bright-to-dark ramp: no particles there.
    img = spheres((80, 120), [(90, 40)], [5], noise=0.01)
    img += np.linspace(0, 1, 120)[None, :]
    mask = np.zeros(img.shape, bool)
    mask[:, :60] = True
    p = ferrospot.detect(
        img, NM_PER_PX, radius_nm=(1.0, 5.0), threshold=0.05, mask=mask
    )
    assert len(p) == 0


def test_build_mask_window_exclude_and_depth():
    img = np.full((100, 100), 0.5)
    img[:, :20] = 0.0  # vacuum on the left
    img[:, 80:] = 1.0  # bright mineral on the right
    mask, depth = ferrospot.build_mask(
        img,
        1.0,
        window=(0.3, 0.7),
        smooth_nm=0,
        erode_nm=5,
        exclude=[(40, 0, 50, 10)],
        vacuum_below=0.2,
        depth_band_nm=(0, 30),
    )
    assert depth[50, 20] == 1 and depth[50, 45] == 26
    cols = np.nonzero(mask[50])[0]
    assert cols.min() == 25  # eroded 5 px from the glass edge at x = 20
    assert cols.max() == 48  # depth (x - 19) below 30 nm
    assert not mask[5, 45] and mask[15, 45]
    with pytest.raises(ValueError, match="vacuum"):
        ferrospot.build_mask(img, 1.0, depth_band_nm=(0, 10))


def test_spacing_stats_square_lattice():
    a, r = 20.0, 3.0  # nm
    p = lattice_particles(10, a, r)
    mask = np.ones((200, 200), bool)  # 1 nm pixels, one lattice cell per particle

    s = ferrospot.spacing_stats(p, mask=mask, section_thickness_nm=100.0)

    assert s.count == 100
    # Only the outermost ring is closer to the image edge (a / 2) than to
    # its neighbour (a).
    assert s.nn_count == 64
    assert s.nn_centre_mean_nm == a and s.nn_gap_median_nm == a - 2 * r
    assert s.overlap_fraction == 0
    assert s.d_mean_nm == 2 * r
    assert np.isclose(s.n_a_per_um2, 1e6 / a**2)
    assert np.isclose(s.clark_evans, 2.0)  # a / (0.5 a)
    n_v = 1 / a**2 / (100 + 2 * r)
    assert np.isclose(s.nn3d_centre_mean_nm, 0.554 * n_v ** (-1 / 3))
    assert np.isclose(s.volume_fraction, n_v * 4 / 3 * np.pi * r**3)


def test_clark_evans_of_poisson_pattern_is_one():
    rng = np.random.default_rng(3)
    n, size = 4000, 1000.0
    x, y = rng.uniform(0, size, (2, n))
    p = ferrospot.Particles(x, y, np.full(n, 0.1), np.ones(n), 1.0)
    s = ferrospot.spacing_stats(p, mask=np.ones((1000, 1000), bool))
    assert s.clark_evans == pytest.approx(1.0, abs=0.04)
    assert s.nn_random_nm == pytest.approx(0.5 / np.sqrt(n / size**2))


def test_border_correction():
    # A pair deep inside, and a pair hugging the left edge of the mask.
    x = np.array([50.0, 56.0, 3.0, 3.0])
    y = np.array([50.0, 50.0, 20.0, 40.0])
    p = ferrospot.Particles(x, y, np.ones(4), np.ones(4), 1.0)
    nb = ferrospot.neighbours(p, np.ones((100, 100), bool))
    assert list(nb["counted"]) == [True, True, False, False]
    assert np.allclose(nb["edge_nm"][2:], 4)  # pixel 3 is 4 pixels from outside
    assert ferrospot.spacing_stats(p, mask=np.ones((100, 100), bool)).nn_count == 2
    assert ferrospot.spacing_stats(p).nn_count == 4  # no mask, no correction


def test_overlap_fraction_and_gaps():
    # Pair one overlaps (gap -1), pair two is 3 nm apart; far from each other.
    x = np.array([0.0, 3.0, 100.0, 107.0])
    p = ferrospot.Particles(x, np.zeros(4), np.full(4, 2.0), np.ones(4), 1.0)
    nb = ferrospot.neighbours(p)
    assert np.allclose(nb["gap_nm"], [-1, -1, 3, 3])
    s = ferrospot.spacing_stats(p)
    assert s.overlap_fraction == 0.5 and s.nn_gap_mean_nm == 1.0


def test_depth_profile():
    # 10 x 10 lattice, 10 nm apart, depth = y; vacuum above y = 0.
    p = lattice_particles(10, 10.0, 1.0)
    depth = np.repeat(np.arange(100.0)[:, None], 100, axis=1)
    rows = ferrospot.depth_profile(p, depth, bin_nm=20, mask=np.ones((100, 100), bool))
    assert [r["depth_lo_nm"] for r in rows] == [0, 20, 40, 60, 80]
    assert [r["count"] for r in rows] == [20] * 5
    assert all(r["area_nm2"] == 2000 for r in rows)
    assert rows[2]["nn_centre_mean_nm"] == 10 and rows[2]["nn_count"] == 16
    assert rows[2]["n_a_per_um2"] == pytest.approx(1e4)


def test_measure_scale_bar():
    img = np.zeros((100, 200))
    img[90:94, 20:120] = 1.0  # 100 px bar
    img[80:86, 150:160] = 1.0  # some "text" next to it
    assert ferrospot.measure_scale_bar(img, 50.0) == 0.5
    assert ferrospot.measure_scale_bar(img, 50.0, region=(0, 85, 200, 100)) == 0.5
    assert ferrospot.measure_scale_bar(1 - img, 50.0, dark=True) == 0.5
    tilted = np.zeros((20, 200))
    tilted[8:10, 20:70] = tilted[9:11, 70:120] = 1.0  # one-pixel step
    assert ferrospot.measure_scale_bar(tilted, 50.0) == 0.5


def test_load_image_roundtrip(tmp_path):
    a = (np.arange(64, dtype=np.uint16).reshape(8, 8) * 1000).astype(np.uint16)
    Image.fromarray(a).save(tmp_path / "x.tif")
    Image.fromarray(np.stack([a // 257] * 3, -1).astype(np.uint8)).save(
        tmp_path / "x.png"
    )
    assert np.allclose(ferrospot.load_image(tmp_path / "x.tif"), a / 65535)
    assert np.allclose(
        ferrospot.load_image(tmp_path / "x.png"), a / 65535, atol=1 / 255
    )


def test_overlay_and_csv(tmp_path):
    img = spheres((60, 60), [(30, 30)], [8])
    p = ferrospot.detect(img, NM_PER_PX, radius_nm=(2.0, 10.0))
    mask = np.zeros(img.shape, bool)
    mask[5:55, 5:55] = True
    rgb = ferrospot.overlay(img, p, mask)
    assert rgb.shape == (60, 60, 3) and (rgb == [255, 0, 0]).all(-1).any()
    assert (rgb[5, 5:55] == [0, 255, 255]).all()
    text = p.to_csv(tmp_path / "p.csv", nn_nm=np.array([1.5]))
    assert text.splitlines()[0] == "x_nm,y_nm,r_nm,contrast,nn_nm"
    assert (tmp_path / "p.csv").read_text() == text


def test_save_note_appends(tmp_path):
    t = datetime(2026, 10, 15, 14, 30, 0, tzinfo=UTC)
    path = ferrospot.save_note(
        tmp_path / "notes",
        "data/grain_A.tif",
        "first look",
        {"threshold": 0.05},
        {"count": 12},
        when=t,
    )
    ferrospot.save_note(
        tmp_path / "notes", "data/grain_A.tif", "second", when=t.replace(minute=45)
    )
    assert path == tmp_path / "notes" / "grain_A_2026-10-15.md"
    text = path.read_text()
    assert text.startswith("# data/grain_A.tif (2026-10-15)")
    assert text.count("# data/grain_A.tif") == 1
    assert "## 14:30:00" in text and "## 14:45:00" in text
    assert "- threshold: 0.05" in text and "- count: 12" in text and "second" in text


def synthetic_config(tmp_path):
    """A 0.5 nm/px image with a burned-in 50 nm bar, vacuum at the top."""
    centres = grid(6, 25, 2, seed=4) + [0, 30]
    img = 0.3 + 0.5 * spheres((220, 200), centres, np.full(36, 5.0), noise=0.02)
    img[:20] = 0.0  # vacuum
    img[205:209, 10:110] = 1.0  # 100 px bar
    Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8)).save(
        tmp_path / "synthetic.png"
    )
    cfg = {
        "name": "synthetic",
        "source": {"path": "synthetic.png", "crop": [0, 0, 200, 200]},
        "scale": {"bar_nm": 50, "box": [0, 200, 200, 220]},
        "scale_note": "burned-in bar",
        "mask": {"vacuum_below": 0.15, "smooth_nm": 0.5},
        "detect": {"detector": "log", "radius_nm": [1.0, 4.0], "threshold": 0.05},
        "section_thickness_nm": 50,
        "depth_bin_nm": 20,
        "expected": {"count": 36, "nn_centre_mean_nm": 12.5},
    }
    ferrospot.save_config(cfg, tmp_path / "synthetic.json")
    return cfg


def test_config_round_trip(tmp_path):
    cfg = synthetic_config(tmp_path)
    assert ferrospot.load_config(tmp_path / "synthetic.json") == cfg
    (tmp_path / "bad.json").write_text(json.dumps({**cfg, "treshold_note": "x"}))
    with pytest.raises(ValueError, match="treshold_note"):
        ferrospot.load_config(tmp_path / "bad.json")


def test_cli_run_and_table(tmp_path, capsys):
    synthetic_config(tmp_path)
    out = tmp_path / "results"
    cli.main(["run", str(tmp_path / "*.json"), "--out", str(out)])

    summary = json.loads((out / "synthetic.json").read_text())
    assert summary["nm_per_px"] == 0.5
    stats = summary["stats"]
    assert stats["count"] == 36
    assert stats["nn_centre_mean_nm"] == pytest.approx(12.5, rel=0.1)
    assert summary["depth_profile"][0]["depth_lo_nm"] == 0
    csv = (out / "synthetic_particles.csv").read_text().splitlines()
    assert csv[0].endswith("nn_nm,gap_nm,edge_nm,counted,depth_nm")
    assert len(csv) == 37
    assert (out / "synthetic_overlay.png").exists()

    capsys.readouterr()
    cli.main(["table", str(out)])
    text = capsys.readouterr().out
    assert "| synthetic | 36 |" in text and "|   expected | 36 |" in text
