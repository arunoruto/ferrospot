"""Masks: where to count particles, and depth below a vacuum interface.

Intensities (``window``, ``vacuum_below``, ``vacuum_above``) are in the
units of the loaded image (0 to 1 for 8- and 16-bit files), before any
bright-field inversion. Boxes are ``(x0, y0, x1, y1)`` in pixels,
end-exclusive.
"""

import numpy as np
from scipy import ndimage as ndi


def build_mask(
    image: np.ndarray,
    nm_per_px: float,
    *,
    roi: tuple[int, int, int, int] | None = None,
    exclude: list | tuple = (),
    window: tuple[float, float] | None = None,
    smooth_nm: float = 1.0,
    erode_nm: float = 0.0,
    vacuum_below: float | None = None,
    vacuum_above: float | None = None,
    depth_band_nm: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Region to count particles in, and the depth of every pixel.

    The mask is the intersection of everything given:

    * ``roi``: inside this box;
    * ``exclude``: outside these boxes (labels, insets, scale bars);
    * ``window``: smoothed intensity within ``[lo, hi]`` (e.g. the glass
      between dark vacuum and a bright mineral), with holes filled (the
      particles in it), then eroded by ``erode_nm`` to stay away from its
      edges;
    * ``vacuum_below`` / ``vacuum_above``: the vacuum is the largest
      connected region of smoothed intensity below (dark-field, HAADF) or
      above (bright field) this level. It is never part of the mask, and
      it defines the depth map: distance to the vacuum in nm;
    * ``depth_band_nm``: depth within ``[d0, d1)``.

    Parameters
    ----------
    image : numpy.ndarray
        2D grayscale image.
    nm_per_px : float
        Pixel size in nm.
    roi, exclude, window, erode_nm, vacuum_below, vacuum_above, depth_band_nm
        See above; all optional.
    smooth_nm : float
        Gaussian sigma applied before the ``window`` and vacuum
        thresholds, so single particles and noise do not punch holes.

    Returns
    -------
    mask : numpy.ndarray of bool
    depth_nm : numpy.ndarray or None
        Distance to the vacuum in nm (0 in the vacuum); None without a
        vacuum level.
    """
    img = np.asarray(image, float)
    mask = np.ones(img.shape, bool)
    smooth = ndi.gaussian_filter(img, smooth_nm / nm_per_px) if smooth_nm else img
    if window is not None:
        lo, hi = window
        # Filling holes keeps bright particles inside the window region.
        mask &= ndi.binary_fill_holes((smooth >= lo) & (smooth <= hi))
        if erode_nm:
            mask &= ndi.distance_transform_edt(mask) * nm_per_px > erode_nm
    if roi is not None:
        mask &= _box(img.shape, roi)
    for box in exclude:
        mask &= ~_box(img.shape, box)

    depth = None
    if vacuum_below is not None or vacuum_above is not None:
        vac = (
            smooth < vacuum_below if vacuum_below is not None else smooth > vacuum_above
        )
        labels, n = ndi.label(vac)
        if n == 0:
            raise ValueError("no vacuum found at that level")
        vac = labels == np.argmax(np.bincount(labels.ravel())[1:]) + 1
        depth = ndi.distance_transform_edt(~vac) * nm_per_px
        mask &= ~vac
        if depth_band_nm is not None:
            d0, d1 = depth_band_nm
            mask &= (depth >= d0) & (depth < d1)
    elif depth_band_nm is not None:
        raise ValueError("depth_band_nm needs vacuum_below or vacuum_above")
    return mask, depth


def _box(shape, box):
    x0, y0, x1, y1 = box
    out = np.zeros(shape, bool)
    out[y0:y1, x0:x1] = True
    return out
