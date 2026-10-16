"""Find iron nanoparticles in TEM / HAADF-STEM images."""

from ferrospot.config import analyse, load_config, save_config
from ferrospot.detect import Particles, detect, normalise
from ferrospot.image import load_image, measure_scale_bar, overlay
from ferrospot.masks import build_mask
from ferrospot.notes import save_note
from ferrospot.stats import SpacingStats, depth_profile, neighbours, spacing_stats

__version__ = "0.1.0"

__all__ = [
    "Particles",
    "SpacingStats",
    "analyse",
    "build_mask",
    "depth_profile",
    "detect",
    "load_config",
    "load_image",
    "measure_scale_bar",
    "neighbours",
    "normalise",
    "overlay",
    "save_config",
    "save_note",
    "spacing_stats",
]
