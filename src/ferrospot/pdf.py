"""Figure images from PDFs, through poppler-utils (``pdfimages``, ``pdftoppm``).

The tools must be on ``PATH``. On NixOS, run inside
``nix shell nixpkgs#poppler-utils`` (the devenv shell has them).

Image numbers are those of ``pdfimages -list -f <page> -l <page>``: they
start at 0 on the page and count soft masks (``smask``) too.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ferrospot.image import load_image


def _run(tool: str, *args) -> str:
    exe = shutil.which(tool)
    if exe is None:
        raise FileNotFoundError(
            f"{tool} not found: install poppler-utils "
            "(NixOS: nix shell nixpkgs#poppler-utils)"
        )
    done = subprocess.run(
        [exe, *map(str, args)], check=True, capture_output=True, text=True
    )
    return done.stdout


def extract_image(pdf: str | Path, page: int, index: int) -> np.ndarray:
    """One embedded image at its native resolution, as grayscale.

    Parameters
    ----------
    pdf : str or Path
        The PDF file.
    page : int
        Page number, from 1.
    index : int
        Image number on that page (see the module docstring).

    Returns
    -------
    numpy.ndarray
        2D float array in [0, 1].
    """
    with tempfile.TemporaryDirectory() as tmp:
        _run("pdfimages", "-png", "-f", page, "-l", page, pdf, f"{tmp}/img")
        path = Path(tmp) / f"img-{index:03d}.png"
        if not path.exists():
            n = len(list(Path(tmp).glob("img-*.png")))
            raise IndexError(f"page {page} of {pdf} has {n} images, no {index}")
        return load_image(path)


def render_page(pdf: str | Path, page: int, dpi: float = 600) -> np.ndarray:
    """A page rendered in grayscale at ``dpi`` (vector scale bars included)."""
    with tempfile.TemporaryDirectory() as tmp:
        _run(
            "pdftoppm", "-r", dpi, "-f", page, "-l", page, "-gray", "-png",
            "-singlefile", pdf, f"{tmp}/page",
        )  # fmt: skip
        return load_image(Path(tmp) / "page.png")


def image_ppi(pdf: str | Path, page: int, index: int) -> float:
    """Horizontal resolution (pixels per inch) an image is placed at.

    A page rendered at ``dpi`` shows the image at ``dpi / ppi`` render
    pixels per image pixel. ``pdfimages -list`` rounds the value to an
    integer, i.e. to about 0.3 % at 300 ppi.
    """
    out = _run("pdfimages", "-list", "-f", page, "-l", page, pdf)
    for line in out.splitlines()[2:]:
        cols = line.split()
        if int(cols[1]) == index:
            return float(cols[12])
    raise IndexError(f"page {page} of {pdf} has no image {index}")
