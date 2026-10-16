"""Particle size by half-maximum region growing.

Each particle's size is the equivalent-circle radius of the connected
region around its centre that is brighter than halfway between its peak
and its local background. This uses only the image, not the detector's
scale, so LoG and DoG detections are sized the same way.

For a sphere whose projected intensity goes with its chord length,
``sqrt(R^2 - rho^2)``, the half-maximum outline has radius
``R sqrt(3) / 2 = 0.87 R``; blur widens it again, so a well-resolved
sphere reads 10 to 15 % small and one only a few pixels across can read
large. The value is reported as measured, without a model correction.
"""

import numpy as np
from scipy import ndimage as ndi

# Smoothing before thresholding, so pixel noise does not fray the
# half-maximum outline: a quarter of the smallest radius, at least 1 px.
SMOOTH_PER_RADIUS = 0.25
# The background is a low percentile of the ring around the particle,
# not its median: between dense particles the ring is mostly particles.
BG_PERCENTILE = 10


def half_max_radii(img, inside, xs, ys, r_guess, r_min, r_max):
    """Half-maximum radius and contrast of each particle.

    Parameters
    ----------
    img : numpy.ndarray
        Rescaled image, particles bright (see :func:`ferrospot.normalise`).
    inside : numpy.ndarray of bool
        Region of interest; the background ring and the half-maximum
        region only use pixels in it.
    xs, ys : numpy.ndarray
        Particle centres in pixels.
    r_guess : numpy.ndarray or float
        Rough radius in pixels: the peak is searched within half of it,
        the background is the ``BG_PERCENTILE`` percentile of the ring
        from 1.5 to 2.5 times it.
    r_min, r_max : float
        Smallest and largest radius in pixels. The image is smoothed with
        a Gaussian of sigma ``SMOOTH_PER_RADIUS * r_min`` (at least 1 px),
        and the region is cut at ``1.5 r_max`` from the centre.

    Returns
    -------
    radius, contrast : numpy.ndarray
        Equivalent-circle radius of the half-maximum region in pixels
        (NaN where the centre is not above the half level), and peak minus
        background in units of ``img``.
    """
    n = len(xs)
    if n == 0:
        return np.empty(0), np.empty(0)
    r_guess = np.broadcast_to(np.asarray(r_guess, float), (n,))
    smooth = ndi.gaussian_filter(img, max(1.0, SMOOTH_PER_RADIUS * r_min))
    iy = np.clip(np.round(ys).astype(int), 0, img.shape[0] - 1)
    ix = np.clip(np.round(xs).astype(int), 0, img.shape[1] - 1)
    # Each pixel belongs to its nearest centre (a Voronoi cell), so the
    # region of one particle cannot leak into a touching neighbour. The
    # cut lies midway, so of two unequal touching particles the larger
    # one reads small. (A watershed cut along the intensity saddle was
    # tried: it loses particles whose centre is not exactly on a maximum.)
    seeds = np.zeros(img.shape, int)
    seeds[iy, ix] = np.arange(1, n + 1)
    near = ndi.distance_transform_edt(
        seeds == 0, return_distances=False, return_indices=True
    )
    owner = seeds[near[0], near[1]]

    radius, contrast = np.full(n, np.nan), np.zeros(n)
    for i in range(n):
        rg = r_guess[i]
        half = max(2.5 * rg, 1.5 * r_max) + 2
        y0, x0 = max(0, int(ys[i] - half)), max(0, int(xs[i] - half))
        win = (slice(y0, int(ys[i] + half) + 1), slice(x0, int(xs[i] + half) + 1))
        s = smooth[win]
        yy, xx = np.indices(s.shape)
        rho = np.hypot(yy + y0 - ys[i], xx + x0 - xs[i])
        core = rho <= max(rg / 2, 1.0)
        ring = (rho >= 1.5 * rg) & (rho <= 2.5 * rg) & inside[win]
        peak_at = np.unravel_index(np.argmax(np.where(core, s, -np.inf)), s.shape)
        peak = s[peak_at]
        bg = np.percentile(s[ring], BG_PERCENTILE) if ring.any() else s.min()
        contrast[i] = peak - bg
        region = (s >= (peak + bg) / 2) & (owner[win] == i + 1) & inside[win]
        region &= rho <= 1.5 * r_max
        labels, _ = ndi.label(region)
        if contrast[i] > 0 and labels[peak_at]:
            radius[i] = np.sqrt((labels == labels[peak_at]).sum() / np.pi)
    return radius, contrast
