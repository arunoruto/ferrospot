"""Synthetic sections: 3D packings of iron spheres projected into TEM images.

A :class:`Section` is a slab of matrix (glass or mineral) of thickness
``t`` holding iron spheres, periodic in x and y; :func:`section` fills one
with a polydisc packing. :func:`simulate_image` projects it along z into
a HAADF / dark-field or bright-field image with the pixel size, blur,
noise and 8-bit quantisation of a published figure, ready for
:func:`ferrospot.detect`.

Image model (projected thickness, Beer-Lambert). Per pixel, the electron
path crosses ``L`` nm of iron (the sum of the sphere chords, clipped to
the section) and ``t - L`` nm of matrix. With ``tau = k_m (t - L) + k_Fe
L``:

* bright field: ``I = exp(-tau)`` (electrons scattered out of the
  objective aperture are lost; particles dark);
* HAADF / dark field: ``S = 1 - exp(-tau)`` (the share scattered onto the
  annular detector; particles bright, saturating for thick, heavy
  columns).

``k`` is an inverse mean free path. At 100 nm and these coefficients the
exponent is small and both are close to linear in mass-thickness, so only
the ratio ``k_Fe / k_m`` matters: about 4 for HAADF (incoherent
scattering grows as ``Z^1.7``; Fe, Z = 26 and 7.9 g/cm^3, against a
silicate of mean Z about 11 and 3 g/cm^3, at similar atom densities) and
about 3 for bright field (mass-thickness contrast). Diffraction
(Bragg) contrast of crystalline iron, why some particles are much darker
than others in bright field, is approximated by an optional random
visibility factor per particle on ``k_Fe`` (``bragg``, see
:data:`BRAGG`). Not modelled: phase / Fresnel contrast, dynamical
scattering and channelling, non-spherical particles, a wedge-shaped
section. The image is then blurred with a Gaussian (probe and figure
resolution), averaged into the panel's pixels, stretched like a
published figure (0.5 to 99.9 percentile), noised and quantised to 8 bit.
"""

from dataclasses import dataclass, replace

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

# Inverse mean free paths in 1/nm: (matrix, iron).
COEFFS = {"haadf": (1e-3, 4e-3), "bf": (1 / 150, 1 / 50)}
# Bragg visibility in bright field: a share p_strong of the particles is
# oriented for strong diffraction (k_Fe times k_strong), p_mid moderately,
# the rest shows mass-thickness contrast only (k_weak). ``haadf`` scales
# the effect in HAADF / dark field (0: none; HAADF is nearly incoherent).
BRAGG = {"p_strong": 0.25, "p_mid": 0.25, "k_strong": 4.0, "k_mid": 2.0,
         "k_weak": 1.0, "haadf": 0.0}  # fmt: skip


@dataclass
class Section:
    """Iron spheres in a slab, periodic in x and y.

    Attributes
    ----------
    x, y, z, r : numpy.ndarray
        Centres and radii in nm; x, y in ``[0, size_nm)``. The section is
        ``0 <= z <= thickness_nm``; spheres reaching out of it are cut, and
        those beyond it are bulk material for the 3D statistics only.
    size_nm : float
        Edge of the periodic square.
    thickness_nm : float
        Section thickness.
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    r: np.ndarray
    size_nm: float
    thickness_nm: float

    def visible(self) -> np.ndarray:
        """Spheres that reach into the section."""
        return (self.z + self.r > 0) & (self.z - self.r < self.thickness_nm)

    def subset(self, keep) -> "Section":
        """The spheres selected by a boolean or index array."""
        return Section(
            self.x[keep], self.y[keep], self.z[keep], self.r[keep],
            self.size_nm, self.thickness_nm,
        )  # fmt: skip


def radius_spec(gmean_nm: float, gstd: float = 1.0, rmin_nm: float = 0.0):
    """Lognormal radii (median ``gmean_nm``), monodisperse for ``gstd == 1``.

    Truncated at two geometric standard deviations (polydisc's default)
    and below at ``rmin_nm``.
    """
    from polydisc import RadiusSpec

    lo, hi = max(rmin_nm, gmean_nm / gstd**2), gmean_nm * gstd**2
    if gstd <= 1.0 or lo >= hi:
        return float(max(gmean_nm, rmin_nm))
    return RadiusSpec.lognormal(gmean_nm, gstd, rmin=lo, rmax=hi)


def section(
    size_nm: float,
    thickness_nm: float,
    radii,
    *,
    gap_nm: float = 0.0,
    keep: float = 1.0,
    packing: str = "gap",
    alpha: float = 1.0,
    depth: dict | None = None,
    seed: int = 0,
) -> Section:
    """Fill a section with a saturated packing, optionally thinned.

    polydisc packs a periodic box, thicker than the section by a margin on
    both sides so its walls do not order the spheres inside, to
    saturation (random sequential addition) with at least ``gap_nm``
    between surfaces; then each sphere is kept with probability ``keep``.
    Saturation at a given gap fixes the density ("regular" states);
    thinning a dense packing gives random hard-core states of any lower
    density. (polydisc's ``max_points`` stops the growth early, which
    leaves a compact cluster, not a sparse uniform packing.)

    ``packing="alpha"`` scales the exclusion with the radii instead:
    centres of spheres i and j at least ``alpha * (r_i + r_j)`` apart
    (``alpha = 1``: touching allowed), i.e. a saturated packing of spheres
    ``alpha`` times larger, shrunk back; ``gap_nm`` is then ignored.

    Parameters
    ----------
    size_nm : float
        Edge of the periodic square (x and y).
    thickness_nm : float
        Section thickness.
    radii : float or polydisc.RadiusSpec
        Radius in nm, or a distribution (see :func:`radius_spec`).
    gap_nm : float
        Minimum surface-to-surface distance.
    keep : float
        Share of spheres kept.
    packing : {"gap", "alpha"}
        Fixed surface gap ``gap_nm``, or radius-scaled exclusion ``alpha``.
    alpha : float
        Exclusion factor for ``packing="alpha"``, at least 1.
    depth : dict, optional
        Number-density profile along z, the grain surface at ``z = 0``
        (no spheres above it): ``{"kind": "exp", "scale_nm": L}`` keeps
        each sphere with probability ``exp(-z / L)``, ``{"kind": "peak",
        "z0_nm": z0, "w_nm": w}`` with ``exp(-(z - z0)^2 / 2 w^2)``.
        Default: uniform.
    seed : int
        Seed of the packing and the thinning.

    Returns
    -------
    Section
    """
    from polydisc import Box, as_radius_spec, sample

    spec, scale = as_radius_spec(radii), 1.0
    if packing == "alpha":
        if alpha < 1:
            raise ValueError(f"alpha must be >= 1, not {alpha}")
        spec, scale, gap_nm = _scaled(spec, alpha), alpha, 0.0
    elif packing != "gap":
        raise ValueError(f"packing must be 'gap' or 'alpha', not {packing!r}")
    margin = 2 * spec.rmax + gap_nm
    box = Box(
        [size_nm, size_nm, thickness_nm + 2 * margin], periodic=(True, True, False)
    )
    pack = sample(box, spec, gap=gap_nm, seed=seed)
    (x, y, z), r = pack.points.T, pack.radii / scale
    sec = Section(x % size_nm, y % size_nm, z - margin, r, size_nm, thickness_nm)
    if depth:
        f = depth_weight(sec.z, depth)
        sec = sec.subset(np.random.default_rng([seed, 2]).random(len(f)) < f)
    if keep < 1.0:
        sec = sec.subset(np.random.default_rng(seed).random(len(r)) < keep)
    return sec


def _scaled(spec, a: float):
    """A polydisc RadiusSpec with every radius multiplied by ``a``."""
    from polydisc.radii import KIND_LOGNORMAL, KIND_NORMAL

    p0, p1 = spec.p0, spec.p1
    if spec.kind == KIND_LOGNORMAL:
        p0 = p0 + np.log(a)
    elif spec.kind == KIND_NORMAL:
        p0, p1 = p0 * a, p1 * a
    return replace(
        spec,
        values=spec.values * a,
        p0=p0,
        p1=p1,
        rmin=spec.rmin * a,
        rmax=spec.rmax * a,
    )


def depth_weight(z: np.ndarray, depth: dict) -> np.ndarray:
    """Relative number density (at most 1) at depth ``z`` below the surface
    ``z = 0``, for a profile as in :func:`section`; 0 above the surface."""
    z = np.asarray(z, float)
    if depth["kind"] == "exp":
        f = np.exp(-z / depth["scale_nm"])
    elif depth["kind"] == "peak":
        f = np.exp(-0.5 * ((z - depth["z0_nm"]) / depth["w_nm"]) ** 2)
    else:
        raise ValueError(f"depth kind must be 'exp' or 'peak', not {depth['kind']!r}")
    return np.where(z < 0, 0.0, f)


def truth(sec: Section) -> dict:
    """True 3D statistics of the spheres whose centre lies in the section.

    Neighbours are searched among all spheres (also those outside the
    section), periodic in x and y.

    Returns
    -------
    dict
        ``n`` (spheres reaching into the section), ``r_mean_nm``,
        ``nn3d_nm`` (mean centre distance to the nearest neighbour),
        ``gap3d_nm`` (mean surface gap), ``gap3d_r`` (``gap3d_nm /
        r_mean_nm``), ``within_1r`` (share of spheres whose gap is below
        their own radius), ``fill`` (volume fraction) and ``n_v_per_um3``.
    """
    t, L = sec.thickness_nm, sec.size_nm
    inside = (sec.z >= 0) & (sec.z <= t)
    if inside.sum() < 2:
        keys = ["r_mean_nm", "nn3d_nm", "gap3d_nm", "gap3d_r", "within_1r", "fill"]
        return {"n": int(sec.visible().sum()), **dict.fromkeys(keys, np.nan),
                "n_v_per_um3": float(inside.sum() / (L * L * t) * 1e9)}  # fmt: skip
    zoff = -sec.z.min() + 1
    pts = np.column_stack([sec.x, sec.y, sec.z + zoff])
    tree = cKDTree(pts, boxsize=[L, L, 10 * (pts[:, 2].max() + 1)])
    dist, idx = tree.query(pts[inside], k=2)
    nn, j = dist[:, 1], idx[:, 1]
    r = sec.r[inside]
    gap = nn - r - sec.r[j]
    vol = L * L * t
    return {
        "n": int(sec.visible().sum()),
        "r_mean_nm": float(r.mean()),
        "nn3d_nm": float(nn.mean()),
        "gap3d_nm": float(gap.mean()),
        "gap3d_r": float(gap.mean() / r.mean()),
        "within_1r": float(np.mean(gap < r)),
        "fill": float(np.sum(4 / 3 * np.pi * r**3) / vol),
        "n_v_per_um3": float(inside.sum() / vol * 1e9),
    }


def iron_path(sec: Section, nm_per_px: float, weights=None) -> np.ndarray:
    """Iron path length (nm) through the section, sampled at pixel centres.

    The image covers the periodic square, ``ceil(size_nm / nm_per_px)``
    pixels on a side; spheres wrap around its edges. ``weights`` (one per
    sphere) multiply each sphere's chord.
    """
    n = int(np.ceil(sec.size_nm / nm_per_px))
    out = np.zeros((n, n))
    t = sec.thickness_nm
    vis = sec.visible()
    w = np.ones(len(sec.r)) if weights is None else np.asarray(weights, float)
    for x, y, z, r, wi in zip(sec.x[vis], sec.y[vis], sec.z[vis], sec.r[vis], w[vis]):
        lo_x, hi_x = (
            int(np.floor((x - r) / nm_per_px)),
            int(np.ceil((x + r) / nm_per_px)),
        )
        lo_y, hi_y = (
            int(np.floor((y - r) / nm_per_px)),
            int(np.ceil((y + r) / nm_per_px)),
        )
        cols, rows = np.arange(lo_x, hi_x + 1), np.arange(lo_y, hi_y + 1)
        dx = (cols + 0.5) * nm_per_px - x
        dy = (rows + 0.5) * nm_per_px - y
        h = np.sqrt(np.clip(r**2 - dy[:, None] ** 2 - dx[None, :] ** 2, 0, None))
        chord = np.clip(np.minimum(z + h, t) - np.maximum(z - h, 0), 0, None)
        out[np.ix_(rows % n, cols % n)] += wi * chord
    return out


def simulate_image(
    sec: Section,
    nm_per_px: float,
    *,
    mode: str = "haadf",
    psf_fwhm_nm: float | None = None,
    noise: float = 0.02,
    oversample: int = 3,
    coeffs: tuple[float, float] | None = None,
    thickness_sd: float = 0.0,
    bragg: dict | bool | None = None,
    seed: int = 0,
) -> np.ndarray:
    """Project a section into a published-looking TEM image.

    Parameters
    ----------
    sec : Section
        The spheres and the section thickness.
    nm_per_px : float
        Pixel size of the image.
    mode : {"haadf", "bf"}
        HAADF / dark field (bright particles) or bright field (dark).
    psf_fwhm_nm : float, optional
        FWHM of the Gaussian blur (probe, resolution, figure resampling).
        Default: 2 pixels.
    noise : float
        Standard deviation of white Gaussian noise, in units of the
        stretched image (0 to 1).
    oversample : int
        Sub-pixels per pixel and axis for the projection and blur.
    coeffs : tuple of float, optional
        Inverse mean free paths ``(k_matrix, k_iron)`` in 1/nm; default
        ``COEFFS[mode]``.
    thickness_sd : float
        Relative standard deviation of the matrix thickness, varying
        smoothly over an eighth of the image (a Gaussian random field):
        the uneven thickness and composition of a real section, which
        sets how much of a figure's grey range the particles get.
    bragg : dict or True, optional
        Diffraction contrast: per sphere a random factor on its ``k_Fe``
        (:func:`bragg_factors`; ``True`` for the defaults :data:`BRAGG`,
        a dict overrides some of them). Default: none.
    seed : int
        Seed of the noise, the thickness field and the Bragg factors.

    Returns
    -------
    numpy.ndarray
        Image in [0, 1], quantised to 256 levels.
    """
    if mode not in COEFFS:
        raise ValueError(f"mode must be 'haadf' or 'bf', not {mode!r}")
    k_m, k_fe = coeffs or COEFFS[mode]
    fine = nm_per_px / oversample
    path = iron_path(sec, fine)
    w = None if bragg is None else bragg_factors(len(sec.r), mode, bragg, seed)
    heavy = path if w is None else iron_path(sec, fine, w)
    rng = np.random.default_rng(seed)
    matrix = np.full(path.shape, float(sec.thickness_nm))
    if thickness_sd:
        field = ndi.gaussian_filter(
            rng.normal(size=path.shape), path.shape[0] / 8, mode="wrap"
        )
        matrix *= 1 + thickness_sd * field / field.std()
    tau = k_m * np.clip(matrix - path, 0, None) + k_fe * heavy
    img = np.exp(-tau) if mode == "bf" else 1.0 - np.exp(-tau)
    fwhm = 2 * nm_per_px if psf_fwhm_nm is None else psf_fwhm_nm
    img = ndi.gaussian_filter(img, fwhm / 2.3548 / fine, mode="wrap")
    n = img.shape[0] // oversample
    img = img[: n * oversample, : n * oversample]
    img = img.reshape(n, oversample, n, oversample).mean(axis=(1, 3))
    lo, hi = np.percentile(img, [0.5, 99.9])
    img = (img - lo) / (hi - lo or 1.0)
    img += rng.normal(0.0, noise, img.shape)
    return np.round(np.clip(img, 0, 1) * 255) / 255


def bragg_factors(n: int, mode: str, bragg=True, seed: int = 0) -> np.ndarray:
    """Per-sphere factors on ``k_Fe`` for diffraction contrast.

    Strong (``k_strong``) with probability ``p_strong``, moderate
    (``k_mid``) with ``p_mid``, else ``k_weak``; in HAADF the excess over
    1 is scaled by ``haadf`` (default 0: no effect). Seeded by ``seed``
    independently of the image noise.
    """
    b = BRAGG | (bragg if isinstance(bragg, dict) else {})
    u = np.random.default_rng([seed, 1]).random(n)
    v = np.where(u < b["p_strong"], b["k_strong"], b["k_weak"])
    v = np.where((u >= b["p_strong"]) & (u < b["p_strong"] + b["p_mid"]), b["k_mid"], v)
    return 1 + (1.0 if mode == "bf" else b["haadf"]) * (v - 1)


def noise_level(image: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Robust standard deviation of the pixel-scale residual.

    The image minus its Gaussian smoothing (sigma 1 px), as 1.4826 times
    the median absolute deviation inside the mask. Measures noise plus
    any structure at the pixel scale; compare like with like (see
    :func:`ferrospot.calibrate.match_noise`).
    """
    res = image - ndi.gaussian_filter(image, 1.0)
    v = res[mask] if mask is not None else res.ravel()
    return float(1.4826 * np.median(np.abs(v - np.median(v))))
