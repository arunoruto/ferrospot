"""Per-image Markdown notes, as written by the dashboard's comment box."""

from datetime import datetime
from pathlib import Path


def save_note(
    notes_dir: str | Path,
    image_name: str,
    comment: str,
    params: dict | None = None,
    stats: dict | None = None,
    when: datetime | None = None,
) -> Path:
    """Append a timestamped entry to ``<notes_dir>/<image-stem>_<date>.md``.

    Parameters
    ----------
    notes_dir : str or Path
        Folder for the notes; created if missing.
    image_name : str
        Image file name; its stem names the notes file.
    comment : str
        Free text.
    params, stats : dict, optional
        Detection parameters and statistics to record with the comment.
    when : datetime, optional
        Time stamp; default now.

    Returns
    -------
    Path
        The notes file.
    """
    when = when or datetime.now().astimezone()
    path = Path(notes_dir) / f"{Path(image_name).stem}_{when:%Y-%m-%d}.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [] if path.exists() else [f"# {image_name} ({when:%Y-%m-%d})", ""]
    lines += [f"## {when:%H:%M:%S}", "", comment.strip(), ""]
    for title, items in (("Parameters", params), ("Statistics", stats)):
        if items:
            lines += [f"**{title}**", ""]
            lines += [f"- {k}: {_fmt(v)}" for k, v in items.items()]
            lines += [""]
    with path.open("a") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _fmt(v):
    return f"{v:.4g}" if isinstance(v, float) else str(v)
