"""Particle detection: LoG blobs or DoG local maxima, sized at half maximum."""

import io
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.feature import blob_log, peak_local_max
from skimage.morphology import h_maxima
from skimage.segmentation import watershed

from ferrospot.sizes import half_max_radii

# HAADF intensity of a sphere goes with its chord length,
# sqrt(R^2 - rho^2). The scale-normalised LoG response of that profile
# peaks at sigma = 0.60 R (a flat disc would peak at R / sqrt(2)).
SIGMA_PER_RADIUS = 0.60
# Number of LoG scales between the smallest and largest radius.
NUM_SCALES = 20
# Ratio of the two Gaussians of the DoG band-pass.
DOG_RATIO = 1.6


@dataclass
class Particles:
    """Detected particles: centres and radii in nm.

    Attributes
    ----------
    x_nm, y_nm : numpy.ndarray
        Centre coordinates in nm, from the top-left image corner
        (x to the right, y down).
    r_nm : numpy.ndarray
        Equivalent-circle radius of the half-maximum region (see
        :mod:`ferrospot.sizes`).
    contrast : numpy.ndarray
        Peak minus local background, in units of the rescaled image
        (0.5 to 99.9 percentile of the region, inverted for bright field).
    nm_per_px : float
        Pixel size the coordinates were converted with.
    """

    x_nm: np.ndarray
    y_nm: np.ndarray
    r_nm: np.ndarray
    contrast: np.ndarray
    nm_per_px: float

    def __len__(self) -> int:
        return len(self.r_nm)

    def subset(self, keep) -> "Particles":
        """The particles selected by a boolean or index array."""
        x, y, r, c = self.x_nm, self.y_nm, self.r_nm, self.contrast
        return Particles(x[keep], y[keep], r[keep], c[keep], self.nm_per_px)

    def to_csv(self, path: str | Path | None = None, **extra) -> str:
        """Write ``x_nm, y_nm, r_nm, contrast`` and any ``extra`` columns.

        Parameters
        ----------
        path : str or Path, optional
            File to write; the text is returned either way.
        **extra : numpy.ndarray
            More per-particle columns, e.g. ``nn_nm=...``.

        Returns
        -------
        str
            The CSV text.
        """
        cols = {f.name: getattr(self, f.name) for f in fields(self)[:4]} | extra
        buf = io.StringIO()
        np.savetxt(
            buf,
            np.column_stack([np.asarray(c, float) for c in cols.values()]),
            delimiter=",",
            header=",".join(cols),
            comments="",
            fmt="%.6g",
        )
        text = buf.getvalue()
        if path is not None:
            Path(path).write_text(text)
        return text


def normalise(image, mask=None, bright_field=False):
    """Rescale an image for detection.

    The 0.5 and 99.9 percentiles inside the mask map to 0 and 1 (clipped),
    a bright-field image is inverted so particles are bright, and every
    pixel outside the mask takes the value of the nearest pixel inside, so
    neither a vacuum edge nor a cut through a thickness ramp looks like a
    particle.

    Parameters
    ----------
    image : numpy.ndarray
        2D grayscale image.
    mask : numpy.ndarray of bool, optional
        Region of interest.
    bright_field : bool
        Particles are dark (bright-field TEM): invert.

    Returns
    -------
    numpy.ndarray
    """
    img = np.asarray(image, dtype=np.float64)
    inside = np.ones(img.shape, bool) if mask is None else np.asarray(mask, bool)
    lo, hi = np.percentile(img[inside], [0.5, 99.9])
    img = np.clip((img - lo) / (hi - lo or 1.0), 0.0, 1.0)
    if bright_field:
        img = 1.0 - img
    if not inside.all():
        near = ndi.distance_transform_edt(
            ~inside, return_distances=False, return_indices=True
        )
        img = img[tuple(near)]
    return img


def detect(
    image: np.ndarray,
    nm_per_px: float,
    *,
    detector: str = "log",
    radius_nm: tuple[float, float] = (1.0, 20.0),
    threshold: float = 0.05,
    mask: np.ndarray | None = None,
    bright_field: bool = False,
    split_touching: bool = False,
) -> Particles:
    """Find roughly round particles and size them at half maximum.

    Parameters
    ----------
    image : numpy.ndarray
        2D grayscale image.
    nm_per_px : float
        Pixel size in nm.
    detector : {"log", "dog"}
        ``"log"``: :func:`skimage.feature.blob_log` over ``NUM_SCALES``
        scales between the two radii; good for mixed sizes.
        ``"dog"``: local maxima of one Difference of Gaussians tuned to
        the smallest radius (sigma ``SIGMA_PER_RADIUS * r_min`` and
        ``DOG_RATIO`` times that), at least ``r_min`` apart; good for
        dense, small particles that LoG merges.
    radius_nm : tuple of float
        Smallest and largest particle radius, in nm. Sizes are grown out
        to at most 1.5 times the largest.
    threshold : float
        Minimum LoG (scale-normalised) or DoG response, in units of the
        rescaled image; lower finds fainter particles and more noise.
    mask : numpy.ndarray of bool, optional
        Region of interest; particles whose centre lies outside are
        dropped.
    bright_field : bool
        Particles are dark on a bright background (bright-field TEM).
    split_touching : bool
        LoG only: look inside every blob for touching particles
        (threshold at half its contrast, watershed on the distance
        transform when it has more than one core).

    Returns
    -------
    Particles
        Only particles whose half-maximum region could be grown.
    """
    shape = np.shape(image)
    inside = np.ones(shape, bool) if mask is None else np.asarray(mask, bool)
    rmin, rmax = (r / nm_per_px for r in radius_nm)
    img = normalise(image, inside, bright_field)

    if detector == "log":
        blobs = blob_log(
            img,
            min_sigma=SIGMA_PER_RADIUS * rmin,
            max_sigma=SIGMA_PER_RADIUS * rmax,
            num_sigma=NUM_SCALES,
            threshold=threshold,
            exclude_border=False,
        )
        y, x, r = blobs[:, 0], blobs[:, 1], blobs[:, 2] / SIGMA_PER_RADIUS
    elif detector == "dog":
        sigma = SIGMA_PER_RADIUS * rmin
        band = ndi.gaussian_filter(img, sigma) - ndi.gaussian_filter(
            img, DOG_RATIO * sigma
        )
        peaks = peak_local_max(
            band,
            min_distance=max(1, round(rmin)),
            threshold_abs=threshold,
            exclude_border=False,
        )
        y, x = peaks.T.astype(float)
        r = np.full(len(y), rmin)
    else:
        raise ValueError(f"detector must be 'log' or 'dog', not {detector!r}")

    keep = inside[y.astype(int), x.astype(int)]
    x, y, r = x[keep], y[keep], r[keep]
    if split_touching and detector == "log":
        x, y, r = _split_all(img, inside, x, y, r, rmin)

    # Sizes are measured after removing background variation on scales
    # well above the particles (thickness gradients, bright-field shading).
    flat = img - ndi.gaussian_filter(img, 2 * rmax)
    radius, contrast = half_max_radii(flat, inside, x, y, r, rmin, rmax)
    ok = np.isfinite(radius)
    s = nm_per_px
    return Particles(x[ok] * s, y[ok] * s, radius[ok] * s, contrast[ok], s)


def _split_all(img, inside, xs, ys, rs, rmin):
    """Replace each LoG blob by its touching parts, if it has several."""
    _, contrast = half_max_radii(img, inside, xs, ys, rs, rmin, rs.max(initial=1.0))
    smooth = ndi.gaussian_filter(img, max(0.5, 0.3 * rmin))
    rows = [
        part
        for x, y, r, c in zip(xs, ys, rs, contrast)
        for part in _split(smooth, x, y, r, smooth[int(y), int(x)] - c / 2, rmin)
    ]
    return np.array(rows, dtype=np.float64).reshape(-1, 3).T


def _split(smooth, x, y, r, level, rmin):
    """Split one blob into touching parts; returns [(x, y, r), ...]."""
    half = 1.2 * r + 1
    y0, x0 = max(0, int(y - half)), max(0, int(x - half))
    win = (slice(y0, int(y + half) + 1), slice(x0, int(x + half) + 1))
    yy, xx = np.indices(smooth[win].shape)
    fg = (smooth[win] > level) & (np.hypot(yy + y0 - y, xx + x0 - x) <= 1.2 * r)
    labels, _ = ndi.label(fg)
    centre = labels[round(y) - y0, round(x) - x0]
    if centre == 0:
        return [(x, y, r)]
    fg = labels == centre
    dist = ndi.distance_transform_edt(fg)
    # h-maxima with h = 1 px ignores the ripples a noisy outline leaves.
    markers, n = ndi.label(h_maxima(dist, 1.0))
    if n < 2:
        return [(x, y, r)]
    regions = watershed(-dist, markers, mask=fg)
    idx = np.arange(1, n + 1)
    # A sphere's half-maximum outline has radius R sqrt(3)/2.
    radii = ndi.maximum(dist, regions, idx) / (np.sqrt(3) / 2)
    cy, cx = np.array(ndi.center_of_mass(fg, regions, idx)).T
    ok = radii >= rmin
    if ok.sum() < 2:
        return [(x, y, r)]
    return list(zip(cx[ok] + x0, cy[ok] + y0, radii[ok]))
