"""Image loading, scale-bar measurement and a quick-look overlay."""

from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.draw import circle_perimeter
from skimage.segmentation import find_boundaries

# Modes PIL can hand over as one numeric channel without conversion.
_NUMERIC_MODES = {"L", "I", "F", "I;16", "I;16B", "I;16L"}


def load_image(path: str | Path) -> np.ndarray:
    """Load a PNG, TIFF or JPEG image as a grayscale float array.

    Parameters
    ----------
    path : str or Path or file-like
        Image file. Colour images are converted to luminance; for a
        multi-frame TIFF only the first frame is read.

    Returns
    -------
    numpy.ndarray
        2D float64 array. Unsigned integer images are scaled to [0, 1]
        by their type's maximum (8- and 16-bit); float and 32-bit integer
        images keep their values.
    """
    with Image.open(path) as im:
        if im.mode not in _NUMERIC_MODES:
            im = im.convert("L")
        a = np.asarray(im)
    if a.dtype.kind == "u":
        return a / np.iinfo(a.dtype).max
    return a.astype(np.float64)


def measure_scale_bar(
    image: np.ndarray,
    bar_nm: float,
    region: tuple[int, int, int, int] | None = None,
    dark: bool = False,
) -> float:
    """Estimate nm per pixel from a horizontal scale bar.

    The image (or ``region``) is thresholded halfway between its minimum
    and maximum, and the widest connected patch of bright (or, with
    ``dark``, dark) pixels is taken as the bar; its horizontal extent is
    its length, so a bar tilted by a degree or two in a scan still
    measures right. This is a
    helper, not a measurement to trust blindly: text next to the bar, end
    ticks, anti-aliasing or a busy background all change the result.
    Check the returned length against the image, and prefer the pixel
    size from the microscope metadata when you have it.

    Parameters
    ----------
    image : numpy.ndarray
        2D grayscale image.
    bar_nm : float
        Length the scale bar stands for, in nm.
    region : tuple of int, optional
        Box ``(x0, y0, x1, y1)`` in pixels (end-exclusive) to search in,
        e.g. the corner holding the bar. Default: the whole image.
    dark : bool
        The bar is dark on a bright background.

    Returns
    -------
    float
        ``bar_nm`` divided by the bar length in pixels.
    """
    if region is not None:
        x0, y0, x1, y1 = region
        image = image[y0:y1, x0:x1]
    bright = image > (image.min() + image.max()) / 2
    if dark:
        bright = ~bright
    if not bright.any() or bright.all():
        raise ValueError("no bar contrast found; is the region right?")
    labels, _ = ndi.label(bright, structure=np.ones((3, 3)))
    widths = [sl[1].stop - sl[1].start for sl in ndi.find_objects(labels)]
    return bar_nm / max(widths)


def overlay(image: np.ndarray, particles, mask: np.ndarray | None = None) -> np.ndarray:
    """Draw the detected particles as red circles, the mask edge in cyan.

    Parameters
    ----------
    image : numpy.ndarray
        2D grayscale image the particles were detected in.
    particles : Particles
        Result of :func:`ferrospot.detect`.
    mask : numpy.ndarray of bool, optional
        Region counted in; its outline is drawn.

    Returns
    -------
    numpy.ndarray
        ``(H, W, 3)`` uint8 RGB image.
    """
    lo, hi = np.percentile(image, [0.5, 99.9])
    gray = (np.clip((image - lo) / (hi - lo or 1), 0, 1) * 255).astype(np.uint8)
    rgb = np.repeat(gray[..., None], 3, axis=-1)
    if mask is not None:
        rgb[find_boundaries(np.asarray(mask, bool), mode="inner")] = (0, 255, 255)
    s = particles.nm_per_px
    for x, y, r in zip(particles.x_nm / s, particles.y_nm / s, particles.r_nm / s):
        rr, cc = circle_perimeter(
            round(y), round(x), max(1, round(r)), shape=image.shape
        )
        rgb[rr, cc] = (255, 0, 0)
    return rgb
