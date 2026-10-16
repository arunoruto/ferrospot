"""Command line: run configs, tabulate results.

``ferrospot run "configs/*.json" --out results`` writes, per config,
``<name>.json`` (summary), ``<name>_particles.csv`` (per particle,
with nearest-neighbour distance, gap, distance to the mask edge, whether
it passed the border correction and, with a vacuum level, depth) and
``<name>_overlay.png``. ``ferrospot table results`` prints a Markdown
table of all summaries, with the ``expected`` values under each row.

``ferrospot calibrate "configs/*.json" --out calibration`` (needs the
``simulate`` extra) simulates each panel's imaging conditions over a grid
of 3D packings (:mod:`ferrospot.calibrate`) and writes, per config,
``<name>_calibration.json`` (reused on the next run), ``<name>_estimate.json``,
``<name>_calibration.png`` and ``<name>_compare.png`` (real crop next to
the simulation at the inferred state), plus ``table.md``. Options:
``--alphas 1,1.25,1.5,2,3`` (radius-scaled hard-core family),
``--bragg`` (bright-field panels: Bragg visibility as a nuisance),
``--depth exp:30`` or ``--depth peak:20:10`` (depth profile, nm) and
``--match nndd`` (fit the NN-distance distribution, not its mean; both
are reported).
"""

import argparse
import glob
import json
from pathlib import Path

from PIL import Image

from ferrospot.config import analyse, load_config, summary
from ferrospot.image import overlay


def run(patterns: list[str], out: str | Path) -> list[Path]:
    """Run every config matching ``patterns``; return the summary files."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    files = sorted({f for p in patterns for f in glob.glob(p)})
    if not files:
        raise FileNotFoundError(f"no config matches {patterns}")
    written = []
    for f in files:
        cfg = load_config(f)
        res = analyse(cfg, Path(f).parent)
        name, p, nb = cfg["name"], res["particles"], res["neighbours"]
        extra = {k: nb[k] for k in ("nn_nm", "gap_nm", "edge_nm", "counted")}
        if res["depth_nm"] is not None:
            iy = (p.y_nm / p.nm_per_px).astype(int)
            ix = (p.x_nm / p.nm_per_px).astype(int)
            extra["depth_nm"] = res["depth_nm"][iy, ix]
        p.to_csv(out / f"{name}_particles.csv", **extra)
        Image.fromarray(overlay(res["image"], p, res["mask"])).save(
            out / f"{name}_overlay.png"
        )
        path = out / f"{name}.json"
        path.write_text(json.dumps(summary(cfg, res), indent=2) + "\n")
        st = res["stats"]
        print(
            f"{name}: n = {st.count}, d = {_pm(st, 'd')}, NN = {_pm(st, 'nn_centre')}"
        )
        written.append(path)
    return written


def calibrate(
    patterns: list[str], out: str | Path, workers: int = 1, *, alphas=(),
    bragg=False, depth: dict | None = None, match: str = "mean",
) -> str:  # fmt: skip
    """2D -> 3D calibration of every config matching ``patterns``.

    A saved calibration is reused if it was made with the same options
    (and, for ``match="nndd"``, recorded the NN-distance quantiles).
    """
    import matplotlib

    matplotlib.use("Agg")
    from ferrospot import calibrate as cb

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    files = sorted({f for p in patterns for f in glob.glob(p)})
    if not files:
        raise FileNotFoundError(f"no config matches {patterns}")
    entries = []
    for f in files:
        cfg = load_config(f)
        name = cfg["name"]
        res = analyse(cfg, Path(f).parent)
        path = out / f"{name}_calibration.json"
        opts = {"alphas": list(alphas), "depth": depth}
        opts["bragg"] = dict(cb.sim.BRAGG) if bragg else None
        cal = cb.Calibration.load(path) if path.exists() else None
        if (
            cal is None
            or cal.options != {k: v for k, v in opts.items() if v}
            or (match == "nndd" and not cal.measured.get("nn_q"))
        ):
            cal = cb.calibrate_panel(cfg, res, workers=workers, **opts)
            cal.save(path)
        est = cb.estimate_3d(res["stats"], cal, match=match)
        other = "nndd" if match == "mean" else "mean"
        if other == "mean" or cal.measured.get("nn_q"):
            est["other"] = cb.estimate_3d(res["stats"], cal, match=other)
        (out / f"{name}_estimate.json").write_text(json.dumps(est, indent=1) + "\n")
        for kind, fig in [
            ("calibration", cb.plot_calibration(cal, est)),
            ("compare", cb.plot_compare(res, cal, est)),
        ]:
            fig.savefig(out / f"{name}_{kind}.png", dpi=110)
            matplotlib.pyplot.close(fig)
        alt = est.get("other", {}).get("gap3d_r")
        print(
            f"{name}: 3D gap {est['gap3d_r']:.1f} r ({match})"
            + (f", {alt:.1f} r ({other})" if alt is not None else ""),
            flush=True,
        )
        entries.append((name, cal, est))
    text = cb.table(entries)
    (out / "table.md").write_text(text + "\n")
    return text


def table(results: str | Path) -> str:
    """Markdown table of every ``<results>/*.json`` summary."""
    head = "| image | n | diameter | NN centre | gap (median) | overlap | Clark-Evans |"
    lines = [head, "|" + "---|" * 7]
    for f in sorted(Path(results).glob("*.json")):
        s = json.loads(f.read_text())
        lines.append(_row(s["config"]["name"], s["stats"]))
        if "expected" in s["config"]:
            lines.append(_row("  expected", s["config"]["expected"]))
    return "\n".join(lines)


def _num(v):
    if v is None:
        return "–"
    return f"{v:.0f}" if abs(v) >= 100 else f"{v:.1f}"


def _pm(st, name):
    get = st.get if isinstance(st, dict) else lambda k: getattr(st, k, None)
    mean, std = get(f"{name}_mean_nm"), get(f"{name}_std_nm")
    return _num(mean) + ("" if std is None else f" ± {_num(std)}")


def _row(name, st):
    ov, ce = st.get("overlap_fraction"), st.get("clark_evans")
    cells = [
        name,
        str(st.get("count", "–")),
        _pm(st, "d"),
        _pm(st, "nn_centre"),
        f"{_pm(st, 'nn_gap')} ({_num(st.get('nn_gap_median_nm'))})",
        "–" if ov is None else f"{100 * ov:.0f} %",
        "–" if ce is None else f"{ce:.2f}",
    ]
    return "| " + " | ".join(cells) + " |"


def _floats(text: str) -> list[float]:
    return [float(v) for v in text.split(",")]


def _depth(text: str) -> dict:
    kind, *v = text.split(":")
    if kind == "exp" and len(v) == 1:
        return {"kind": "exp", "scale_nm": float(v[0])}
    if kind == "peak" and len(v) == 2:
        return {"kind": "peak", "z0_nm": float(v[0]), "w_nm": float(v[1])}
    raise argparse.ArgumentTypeError("use exp:<scale_nm> or peak:<z0_nm>:<w_nm>")


def main(argv: list[str] | None = None) -> None:
    """Entry point of ``ferrospot``."""
    parser = argparse.ArgumentParser(
        prog="ferrospot", description=__doc__.split("\n")[0]
    )
    sub = parser.add_subparsers(dest="command", required=True)
    p_run = sub.add_parser("run", help="run configs")
    p_run.add_argument("configs", nargs="+", help="config files or glob patterns")
    p_run.add_argument("--out", default="results", help="output folder")
    p_table = sub.add_parser("table", help="tabulate results")
    p_table.add_argument("results", nargs="?", default="results")
    p_cal = sub.add_parser("calibrate", help="2D -> 3D calibration by simulation")
    p_cal.add_argument("configs", nargs="+", help="config files or glob patterns")
    p_cal.add_argument("--out", default="calibration", help="output folder")
    p_cal.add_argument("--workers", type=int, default=1, help="processes")
    p_cal.add_argument(
        "--alphas", type=_floats, default=(),
        help="radius-scaled hard-core family, e.g. 1,1.25,1.5,2,3",
    )  # fmt: skip
    p_cal.add_argument(
        "--bragg", action="store_true", help="bright field: Bragg visibility nuisance"
    )
    p_cal.add_argument(
        "--depth", type=_depth, default=None,
        help="depth profile: exp:<scale_nm> or peak:<z0_nm>:<w_nm>",
    )  # fmt: skip
    p_cal.add_argument("--match", choices=("mean", "nndd"), default="mean")
    args = parser.parse_args(argv)
    if args.command == "run":
        run(args.configs, args.out)
    elif args.command == "calibrate":
        print(
            calibrate(
                args.configs, args.out, args.workers, alphas=args.alphas,
                bragg=args.bragg, depth=args.depth, match=args.match,
            )
        )  # fmt: skip
    else:
        print(table(args.results))


if __name__ == "__main__":
    main()
