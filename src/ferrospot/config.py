"""One JSON config per analysed image, and the analysis it describes.

A config records where the image comes from, how its scale was found,
the mask and the detection settings, with free-text ``*_note`` fields for
every choice that is not obvious::

    {
      "name": "example_fig2",
      "description": "Example paper, Fig. 2: nanoparticles in a grain rim",
      "source": {"pdf": "~/papers/example.pdf", "page": 3, "image": 6,
                 "crop": [0, 0, 514, 513]},
      "scale": {"bar_nm": 20, "bar_px": 68},
      "scale_note": "20 nm vector bar, 68 px in a 300 dpi page render",
      "mask": {"vacuum_below": 0.1, "smooth_nm": 2},
      "mask_note": "...",
      "detect": {"detector": "dog", "radius_nm": [0.6, 3], "threshold": 0.02},
      "threshold_note": "...",
      "section_thickness_nm": 100,
      "depth_bin_nm": 20,
      "expected": {"count": 700, "d_mean_nm": 2.5},
      "expected_note": "..."
    }

``source``: either ``path`` (an image file) or ``pdf`` + ``page`` +
``image`` (see :mod:`ferrospot.pdf`); optional ``crop`` box ``(x0, y0,
x1, y1)``. Relative paths are relative to the config file; ``~`` is
expanded.

``scale``, one of:

* ``nm_per_px``;
* ``bar_nm`` + ``bar_px``: a bar measured by hand;
* ``bar_nm`` + ``box``: a bar burned into the image, measured with
  :func:`ferrospot.measure_scale_bar` in this box of the uncropped image
  (``"dark": true`` for a dark bar);
* ``bar_nm`` + ``page_box`` (+ ``dpi``, default 600): a vector bar drawn
  on the page, measured in this box of a page render and converted with
  the image's placement resolution (:func:`ferrospot.pdf.image_ppi`).

``mask`` and ``detect`` take the keyword arguments of
:func:`ferrospot.build_mask` and :func:`ferrospot.detect`; all mask
coordinates refer to the cropped image. ``expected`` holds reference
values (fields of :class:`ferrospot.SpacingStats`) to compare against.
"""

import inspect
import json
import re
from dataclasses import asdict
from pathlib import Path

import numpy as np

from ferrospot import pdf
from ferrospot.detect import detect
from ferrospot.image import load_image, measure_scale_bar
from ferrospot.masks import build_mask
from ferrospot.stats import depth_profile, neighbours, spacing_stats


def _kwargs(func):
    return {
        p.name
        for p in inspect.signature(func).parameters.values()
        if p.kind is p.KEYWORD_ONLY
    }


KEYS = {
    "": {
        "name", "description", "source", "scale", "scale_note", "mask",
        "mask_note", "detect", "threshold_note", "section_thickness_nm",
        "depth_bin_nm", "expected", "expected_note",
    },
    "source": {"path", "pdf", "page", "image", "crop"},
    "scale": {"nm_per_px", "bar_nm", "bar_px", "box", "page_box", "dpi", "dark"},
    "mask": _kwargs(build_mask),
    "detect": _kwargs(detect),
}  # fmt: skip


def load_config(path: str | Path) -> dict:
    """Read a config and check its keys; ``name`` defaults to the file stem."""
    cfg = json.loads(Path(path).read_text())
    for section, allowed in KEYS.items():
        unknown = set(cfg.get(section, {}) if section else cfg) - allowed
        if unknown:
            raise ValueError(f"{path}: unknown {section or 'top-level'} keys {unknown}")
    cfg.setdefault("name", Path(path).stem)
    return cfg


def save_config(cfg: dict, path: str | Path) -> None:
    """Write a config as indented JSON, short lists (boxes) on one line."""
    text = json.dumps(cfg, indent=2, ensure_ascii=False)
    text = re.sub(r"\[\s+([^\[\]{}]*?)\s+\]", _one_line, text)
    Path(path).write_text(text + "\n")


def _one_line(m):
    return "[" + ", ".join(v.strip() for v in m.group(1).split(",")) + "]"


def _path(p, base):
    p = Path(p).expanduser()
    return p if p.is_absolute() else Path(base) / p


def nm_per_px(cfg: dict, image: np.ndarray, base: str | Path = ".") -> float:
    """Pixel size from the ``scale`` section (``image`` uncropped)."""
    sc, src = cfg["scale"], cfg["source"]
    if "nm_per_px" in sc:
        return float(sc["nm_per_px"])
    if "bar_px" in sc:
        return sc["bar_nm"] / sc["bar_px"]
    if "box" in sc:
        return measure_scale_bar(image, sc["bar_nm"], sc["box"], sc.get("dark", False))
    if "page_box" in sc:
        file, dpi = _path(src["pdf"], base), sc.get("dpi", 600)
        page = pdf.render_page(file, src["page"], dpi)
        per_render_px = measure_scale_bar(
            page, sc["bar_nm"], sc["page_box"], sc.get("dark", False)
        )
        return per_render_px * dpi / pdf.image_ppi(file, src["page"], src["image"])
    raise ValueError("scale needs nm_per_px, bar_px, box or page_box")


def analyse(cfg: dict, base: str | Path = ".") -> dict:
    """Run one config: load, scale, crop, mask, detect, measure.

    Parameters
    ----------
    cfg : dict
        A config (see the module docstring).
    base : str or Path
        Folder that relative source paths are relative to.

    Returns
    -------
    dict
        ``image`` (cropped), ``nm_per_px``, ``mask``, ``depth_nm`` (or
        None), ``particles``, ``neighbours`` (per-particle arrays from
        :func:`ferrospot.stats.neighbours`), ``stats``
        (:class:`ferrospot.SpacingStats`) and ``depth_profile`` (list, empty
        without a vacuum level).
    """
    src = cfg["source"]
    if "path" in src:
        image = load_image(_path(src["path"], base))
    else:
        image = pdf.extract_image(_path(src["pdf"], base), src["page"], src["image"])
    s = nm_per_px(cfg, image, base)
    if "crop" in src:
        x0, y0, x1, y1 = src["crop"]
        image = image[y0:y1, x0:x1]
    mask, depth = build_mask(image, s, **cfg.get("mask", {}))
    p = detect(image, s, mask=mask, **cfg.get("detect", {}))
    t = cfg.get("section_thickness_nm")
    profile = []
    if depth is not None:
        profile = depth_profile(p, depth, cfg.get("depth_bin_nm", 20.0), mask)
    return {
        "image": image,
        "nm_per_px": s,
        "mask": mask,
        "depth_nm": depth,
        "particles": p,
        "neighbours": neighbours(p, mask),
        "stats": spacing_stats(p, mask=mask, section_thickness_nm=t),
        "depth_profile": profile,
    }


def summary(cfg: dict, result: dict) -> dict:
    """JSON-ready summary: the config, pixel size, stats and depth profile."""
    return {
        "config": cfg,
        "nm_per_px": result["nm_per_px"],
        "stats": asdict(result["stats"]),
        "depth_profile": result["depth_profile"],
    }
