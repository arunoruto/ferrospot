"""Size and spacing statistics of detected particles.

Assumptions, read before quoting a number:

* The image is a 2D projection of a section of thickness ``t``. Particles
  at different depths overlap in projection, so 2D nearest-neighbour
  distances are *shorter* than the true 3D ones, and 2D gaps can be
  negative ("overlapping in projection").
* Detection bias works the other way: merged neighbours and missed small
  or faint particles make particles look larger and farther apart.
* Border correction: a particle's nearest neighbour only counts if the
  particle is farther from the mask edge (or image edge) than from that
  neighbour; otherwise a closer neighbour might sit just outside.
  Without a mask there is no correction.
* Clark-Evans ratio: the observed mean nearest-neighbour distance over
  ``0.5 / sqrt(N_A)``, the mean for a Poisson (random) pattern of the
  same density. Above 1 is more regular than random (e.g. hard-sphere
  exclusion), below 1 clustered.
* The 3D estimate assumes particles placed independently at random
  (Poisson) in the section. A particle shows up if its centre lies within
  the section or less than one radius outside it, so the number per area
  is ``N_A = n_V (t + d_mean)``. The mean 3D nearest-neighbour centre
  distance of a Poisson process is ``Gamma(4/3) (4 pi n_V / 3)^(-1/3)
  = 0.554 n_V^(-1/3)``; it ignores that particles cannot overlap, so it
  is only a rough guide once the volume fraction reaches a few percent.
"""

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from ferrospot.detect import Particles


@dataclass
class SpacingStats:
    """Summary statistics; ``None`` where the inputs do not allow it.

    Distances are in nm. ``d`` is the half-maximum diameter, ``nn_centre``
    the distance to the nearest neighbour's centre, ``nn_gap`` the
    surface-to-surface gap ``nn_centre - r_i - r_j``; ``nn_count``
    particles pass the border correction and enter the ``nn_*``,
    ``overlap_fraction`` (share with ``nn_gap < 0``) and ``clark_evans``
    fields. ``n_a`` is per square micrometre, ``n_v`` per cubic
    micrometre. The ``nn3d_*``, ``n_v`` and ``volume_fraction`` fields are
    the random-placement estimate described in the module docstring.
    """

    count: int
    area_nm2: float | None
    n_a_per_um2: float | None
    d_mean_nm: float | None
    d_std_nm: float | None
    d_p10_nm: float | None
    d_median_nm: float | None
    d_p90_nm: float | None
    nn_count: int
    nn_centre_mean_nm: float | None
    nn_centre_std_nm: float | None
    nn_centre_p10_nm: float | None
    nn_centre_median_nm: float | None
    nn_centre_p90_nm: float | None
    nn_gap_mean_nm: float | None
    nn_gap_std_nm: float | None
    nn_gap_p10_nm: float | None
    nn_gap_median_nm: float | None
    nn_gap_p90_nm: float | None
    overlap_fraction: float | None
    nn_random_nm: float | None
    clark_evans: float | None
    section_thickness_nm: float | None
    n_v_per_um3: float | None
    nn3d_centre_mean_nm: float | None
    nn3d_gap_mean_nm: float | None
    volume_fraction: float | None


def neighbours(
    particles: Particles,
    mask: np.ndarray | None = None,
    nm_per_px: float | None = None,
) -> dict[str, np.ndarray]:
    """Nearest neighbour of every particle, with the border correction.

    Parameters
    ----------
    particles : Particles
        Particles to pair up (all of them; filter to the mask first).
    mask : numpy.ndarray of bool, optional
        Region counted in; enables the border correction.
    nm_per_px : float, optional
        Pixel size of ``mask``. Default: ``particles.nm_per_px``.

    Returns
    -------
    dict of numpy.ndarray
        ``nn_nm`` (centre distance), ``gap_nm``, ``edge_nm`` (distance to
        the mask or image edge; inf without a mask) and ``counted``
        (``edge_nm > nn_nm``), one value per particle; ``nn_nm`` and
        ``gap_nm`` are NaN for fewer than two particles.
    """
    n = len(particles)
    nn, gap = np.full(n, np.nan), np.full(n, np.nan)
    if n >= 2:
        xy = np.column_stack([particles.x_nm, particles.y_nm])
        dist, idx = cKDTree(xy).query(xy, k=2)
        nn = dist[:, 1]
        gap = nn - particles.r_nm - particles.r_nm[idx[:, 1]]
    edge = np.full(n, np.inf)
    if mask is not None:
        s = nm_per_px or particles.nm_per_px
        # Padding makes the image border an edge, too.
        dt = ndi.distance_transform_edt(np.pad(np.asarray(mask, bool), 1))
        dt = dt[1:-1, 1:-1]
        edge = dt[_pixels(particles, dt.shape, s)] * s
    return {"nn_nm": nn, "gap_nm": gap, "edge_nm": edge, "counted": edge > nn}


def spacing_stats(
    particles: Particles,
    mask: np.ndarray | None = None,
    nm_per_px: float | None = None,
    section_thickness_nm: float | None = None,
) -> SpacingStats:
    """Count, size distribution and nearest-neighbour spacing.

    Parameters
    ----------
    particles : Particles
        Result of :func:`ferrospot.detect`.
    mask : numpy.ndarray of bool, optional
        Region the particles were counted in. Needed for the areal
        density, the border correction, the Clark-Evans ratio and the 3D
        estimate; particles outside it are ignored.
    nm_per_px : float, optional
        Pixel size of ``mask``. Default: ``particles.nm_per_px``.
    section_thickness_nm : float, optional
        Section (lamella) thickness; enables the 3D estimate.

    Returns
    -------
    SpacingStats
    """
    s = nm_per_px or particles.nm_per_px
    area = None
    if mask is not None:
        mask = np.asarray(mask, bool)
        particles = particles.subset(mask[_pixels(particles, mask.shape, s)])
        area = float(mask.sum()) * s**2
    n = len(particles)
    d = 2 * particles.r_nm
    nb = neighbours(particles, mask, s)
    ok = nb["counted"]
    out = dict.fromkeys(SpacingStats.__dataclass_fields__)
    out.update(
        count=n,
        area_nm2=area,
        section_thickness_nm=section_thickness_nm,
        nn_count=int(ok.sum()),
        **_describe("d", d),
        **_describe("nn_centre", nb["nn_nm"][ok]),
        **_describe("nn_gap", nb["gap_nm"][ok]),
    )
    if ok.any():
        out.update(overlap_fraction=np.mean(nb["gap_nm"][ok] < 0))
    if area and n:
        n_a = n / area  # per nm^2
        out.update(n_a_per_um2=n_a * 1e6, nn_random_nm=0.5 / np.sqrt(n_a))
        if ok.any():
            out.update(clark_evans=nb["nn_nm"][ok].mean() / out["nn_random_nm"])
        if section_thickness_nm is not None:
            n_v = n_a / (section_thickness_nm + d.mean())  # per nm^3
            nn3d = 0.554 * n_v ** (-1 / 3)
            out.update(
                n_v_per_um3=n_v * 1e9,
                nn3d_centre_mean_nm=nn3d,
                nn3d_gap_mean_nm=nn3d - d.mean(),
                volume_fraction=n_v * np.mean(4 / 3 * np.pi * particles.r_nm**3),
            )
    return SpacingStats(**{k: _plain(v) for k, v in out.items()})


def depth_profile(
    particles: Particles,
    depth_nm: np.ndarray,
    bin_nm: float = 20.0,
    mask: np.ndarray | None = None,
    nm_per_px: float | None = None,
) -> list[dict]:
    """Statistics in bins of depth below the vacuum interface.

    Nearest neighbours are found among all particles in the mask (with the
    border correction), then each particle is binned by the depth of its
    centre.

    Parameters
    ----------
    particles : Particles
        Result of :func:`ferrospot.detect`.
    depth_nm : numpy.ndarray
        Depth of every pixel in nm, e.g. from :func:`ferrospot.build_mask`.
    bin_nm : float
        Bin width in nm.
    mask : numpy.ndarray of bool, optional
        Region counted in. Default: everywhere.
    nm_per_px : float, optional
        Pixel size of the maps. Default: ``particles.nm_per_px``.

    Returns
    -------
    list of dict
        One row per bin up to the deepest masked pixel: ``depth_lo_nm``,
        ``depth_hi_nm``, ``count``, ``area_nm2``, ``n_a_per_um2``,
        ``d_mean_nm``, ``nn_count``, ``nn_centre_mean_nm``,
        ``nn_gap_mean_nm``, ``overlap_fraction`` (None where empty).
    """
    s = nm_per_px or particles.nm_per_px
    mask = np.ones(depth_nm.shape, bool) if mask is None else np.asarray(mask, bool)
    particles = particles.subset(mask[_pixels(particles, mask.shape, s)])
    nb = neighbours(particles, mask, s)
    depth = depth_nm[_pixels(particles, mask.shape, s)]
    rows = []
    for lo in np.arange(0.0, depth_nm[mask].max(initial=0.0), bin_nm):
        hi = lo + bin_nm
        here = (depth >= lo) & (depth < hi)
        ok = here & nb["counted"]
        area = float(((depth_nm >= lo) & (depth_nm < hi) & mask).sum()) * s**2
        gaps = nb["gap_nm"][ok]
        row = {
            "depth_lo_nm": lo,
            "depth_hi_nm": hi,
            "count": int(here.sum()),
            "area_nm2": area,
            "n_a_per_um2": here.sum() / area * 1e6 if area else None,
            "d_mean_nm": 2 * particles.r_nm[here].mean() if here.any() else None,
            "nn_count": int(ok.sum()),
            "nn_centre_mean_nm": nb["nn_nm"][ok].mean() if ok.any() else None,
            "nn_gap_mean_nm": gaps.mean() if ok.any() else None,
            "overlap_fraction": np.mean(gaps < 0) if ok.any() else None,
        }
        rows.append({k: _plain(v) for k, v in row.items()})
    return rows


def _pixels(particles, shape, s):
    """Pixel indices (rows, cols) of the particle centres."""
    iy = np.clip((particles.y_nm / s).astype(int), 0, shape[0] - 1)
    ix = np.clip((particles.x_nm / s).astype(int), 0, shape[1] - 1)
    return iy, ix


def _describe(name, v):
    """Mean, std, 10th / 50th / 90th percentile of ``v`` (None if empty)."""
    keys = [f"{name}_{k}_nm" for k in ("mean", "std", "p10", "median", "p90")]
    if not len(v):
        return dict.fromkeys(keys)
    return dict(zip(keys, [v.mean(), v.std(), *np.percentile(v, [10, 50, 90])]))


def _plain(v):
    """NumPy scalars to Python numbers, so the stats print and serialise."""
    return v.item() if isinstance(v, np.generic) else v
