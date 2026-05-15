#!/usr/bin/env python3
"""Generate a self-contained PDF validation report from a CDSS Monte Carlo campaign.

Reads results.jsonl (and optionally manifest.json) produced by cdss_monte_carlo.py
and emits a multi-page A4 PDF containing BER/FER/CRC curves, timing analysis, and a
summary table suitable for scientific distribution alongside the binary release.

Usage:
    python3 generate_report.py \\
        --input  /path/to/results.jsonl \\
        --output cdss-validation-report.pdf \\
        [--manifest /path/to/manifest.json] \\
        [--version v1.0.0] \\
        [--title  "CDSS Modem Validation Report"]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import textwrap
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ── Layout / style ────────────────────────────────────────────────────────────

A4 = (8.27, 11.69)

PROFILE_ORDER = ["aligned_awgn", "sync_awgn", "cfo_drift", "tone", "clock", "mixed"]

PROFILE_LABELS: dict[str, str] = {
    "aligned_awgn": "Aligned AWGN",
    "sync_awgn": "Sync AWGN",
    "cfo_drift": "CFO + Drift",
    "tone": "Tone Interference",
    "clock": "Clock Skew",
    "mixed": "Mixed Impairments",
}

_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
_MARKERS = ["o", "s", "^", "D", "v", "P"]

# ── Scalar helpers ────────────────────────────────────────────────────────────


def _f(v: Any) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def _i(v: Any, default: int = 0) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        try:
            return int(float(v))
        except (TypeError, ValueError):
            return default


def _b(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return v.strip().lower() in {"1", "true", "yes", "y", "on"}
    return bool(v) if v is not None else False


def _ok(v: Any) -> bool:
    return isinstance(v, (int, float)) and math.isfinite(float(v))


def _mean(vals: list[float]) -> float:
    fs = [v for v in vals if _ok(v)]
    return sum(fs) / len(fs) if fs else math.nan


def _pct(vals: list[float], p: float) -> float:
    fs = sorted(v for v in vals if _ok(v))
    if not fs:
        return math.nan
    rank = (len(fs) - 1) * p / 100.0
    lo, hi = int(rank), min(int(rank) + 1, len(fs) - 1)
    return fs[lo] + (rank - lo) * (fs[hi] - fs[lo])


def _wilson(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return math.nan, math.nan
    p = successes / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    m = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, c - m), min(1.0, c + m)


# ── Record normalisation ──────────────────────────────────────────────────────


def _norm(raw: dict[str, Any], row: int) -> dict[str, Any]:
    exp = raw.get("experiment") or {}
    sig = raw.get("signal") or {}
    ch = raw.get("channel") or {}
    acq = raw.get("acquisition") or {}
    dec = raw.get("decoder") or {}
    pay = raw.get("payload") or {}
    perf = raw.get("performance") or {}

    profile = str(
        exp.get("profile") or ch.get("profile") or raw.get("profile") or "unknown"
    )
    snr_db = _f(ch.get("snr_db") or raw.get("snr_db"))
    payload_bits = _i(
        sig.get("payload_bits") or raw.get("payload_bits") or pay.get("total_bits", -1),
        -1,
    )
    exit_code = _i(exp.get("exit_code") or raw.get("exit_code", 0))

    crc_ok = _b(dec.get("crc_ok"))
    decode_attempted = _b(dec.get("decode_attempted"))
    sync_search = _b(acq.get("sync_search"))
    frame_detected = _b(acq.get("frame_detected"))
    sync_ok = _b(acq.get("sync_ok") or (frame_detected if sync_search else True))

    bit_errors = _f(pay.get("bit_errors") or raw.get("bit_errors"))
    total_bits = _f(
        pay.get("total_bits") or (payload_bits if payload_bits > 0 else math.nan)
    )
    ber = _f(pay.get("ber") or pay.get("ber_frame") or raw.get("ber"))
    if not _ok(ber) and _ok(bit_errors) and _ok(total_bits) and total_bits > 0:
        ber = bit_errors / total_bits

    iters = _f(
        dec.get("refinement_iters_used")
        or dec.get("iterations_used")
        or dec.get("ldpc_iterations")
    )

    return {
        "row": row,
        "profile": profile,
        "snr_db": snr_db,
        "payload_bits": payload_bits,
        "exit_code": exit_code,
        "crc_ok": crc_ok,
        "decode_attempted": decode_attempted,
        "sync_search": sync_search,
        "sync_ok": sync_ok,
        "ber": ber,
        "bit_errors": bit_errors,
        "total_bits": total_bits,
        "runtime_ms": _f(exp.get("runtime_ms") or perf.get("runtime_ms")),
        "timing_err": _f(acq.get("timing_error_samples")),
        "acq_score": _f(acq.get("acquisition_score")),
        "iterations": iters,
    }


def load_records(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for row, line in enumerate(fh, 1):
            text = line.strip()
            if not text:
                continue
            try:
                records.append(_norm(json.loads(text), row))
            except json.JSONDecodeError:
                pass
    return records


# ── Statistics ────────────────────────────────────────────────────────────────


def _stats(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    if n == 0:
        return {"n": 0}

    crc_ok = sum(1 for r in records if r["crc_ok"])
    cmd_fail = sum(1 for r in records if r["exit_code"] != 0)
    t_errs = sum(r["bit_errors"] for r in records if _ok(r["bit_errors"]))
    t_bits = sum(r["total_bits"] for r in records if _ok(r["total_bits"]))

    ber_frames = [r["ber"] for r in records if _ok(r["ber"])]
    runtimes = [r["runtime_ms"] for r in records if _ok(r["runtime_ms"])]

    sync_recs = [r for r in records if r["sync_search"]]
    sync_ok = sum(1 for r in sync_recs if r["sync_ok"])
    timing_abs = [abs(r["timing_err"]) for r in sync_recs if _ok(r["timing_err"])]
    acq_scores = [r["acq_score"] for r in sync_recs if _ok(r["acq_score"])]

    crc_lo, crc_hi = _wilson(crc_ok, n)
    ber_lo, ber_hi = _wilson(int(t_errs), int(t_bits))

    return {
        "n": n,
        "cmd_fail": cmd_fail,
        "crc_ok": crc_ok,
        "crc_rate": crc_ok / n,
        "crc_ci95_lo": crc_lo,
        "crc_ci95_hi": crc_hi,
        "fer": 1.0 - crc_ok / n,
        "ber_agg": t_errs / t_bits if t_bits > 0 else math.nan,
        "ber_ci95_lo": ber_lo,
        "ber_ci95_hi": ber_hi,
        "ber_mean": _mean(ber_frames),
        "ber_p95": _pct(ber_frames, 95),
        "sync_frames": len(sync_recs),
        "sync_rate": sync_ok / len(sync_recs) if sync_recs else math.nan,
        "timing_mean": _mean(timing_abs),
        "timing_p95": _pct(timing_abs, 95),
        "acq_score_mean": _mean(acq_scores),
        "acq_score_p05": _pct(acq_scores, 5),
        "runtime_p50": _pct(runtimes, 50),
        "runtime_p95": _pct(runtimes, 95),
    }


def _group(
    records: list[dict[str, Any]],
) -> dict[str, dict[int, dict[float, list[dict]]]]:
    """profile → payload_bits → snr_db → [records]"""
    out: dict[str, dict[int, dict[float, list[dict]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for r in records:
        if _ok(r["snr_db"]) and r["payload_bits"] > 0:
            out[r["profile"]][r["payload_bits"]][r["snr_db"]].append(r)
    return out


# ── Figure builders ───────────────────────────────────────────────────────────


def _grid(n: int) -> tuple[plt.Figure, list[plt.Axes]]:
    cols = 2
    rows = math.ceil(n / cols)
    fig, ax_arr = plt.subplots(rows, cols, figsize=A4)
    flat: list[plt.Axes] = list(ax_arr.flat) if hasattr(ax_arr, "flat") else [ax_arr]
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig, flat[:n]


def _page_title(manifest: dict[str, Any], version: str, title: str) -> plt.Figure:
    fig = plt.figure(figsize=A4)
    fig.patch.set_facecolor("white")

    cfg = manifest
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    snrs = cfg.get("snrs", [])
    pays = cfg.get("payloads", [])
    profs = cfg.get("profiles", [])

    meta = [
        ("Version", version or "—"),
        ("Report date", now),
        ("Campaign preset", cfg.get("preset", "—")),
        ("Trials / point", str(cfg.get("trials_per_point", "—"))),
        ("Workers", str(cfg.get("workers", "—"))),
        ("Random seed", str(cfg.get("seed_base", "—"))),
        ("Campaign start", str(manifest.get("created_utc", "—"))),
        ("Payload sizes", ", ".join(f"{p} bits" for p in pays) if pays else "—"),
        (
            "SNR range",
            f"{min(snrs):.0f} … {max(snrs):.0f} dB  ({len(snrs)} points)"
            if snrs
            else "—",
        ),
        ("Channel profiles", ", ".join(profs) if profs else "—"),
    ]

    fig.text(
        0.50,
        0.95,
        title,
        ha="center",
        va="top",
        fontsize=18,
        fontweight="bold",
        color="#1a1a2e",
    )
    fig.text(
        0.50,
        0.91,
        "CHIRP-DSSS-Polar Digital Modem — Monte Carlo Validation",
        ha="center",
        va="top",
        fontsize=11,
        color="#555555",
    )

    fig.add_artist(
        plt.Line2D(
            [0.08, 0.92],
            [0.89, 0.89],
            transform=fig.transFigure,
            color="#cccccc",
            linewidth=0.8,
        )
    )

    y = 0.84
    value_wrap = 50
    for label, value in meta:
        wrapped = textwrap.wrap(str(value), width=value_wrap) or ["—"]
        fig.text(
            0.12,
            y,
            f"{label}:",
            ha="left",
            va="top",
            fontsize=10,
            fontweight="bold",
            color="#333333",
        )
        fig.text(
            0.40,
            y,
            "\n".join(wrapped),
            ha="left",
            va="top",
            fontsize=10,
            color="#111111",
            linespacing=1.3,
        )
        y -= 0.047 * max(1, len(wrapped))

    fig.add_artist(
        plt.Line2D(
            [0.08, 0.92],
            [y - 0.01, y - 0.01],
            transform=fig.transFigure,
            color="#cccccc",
            linewidth=0.8,
        )
    )

    footnote = (
        "Validation methodology: each trial invokes the compiled cdss binary via "
        "`cdss simulate`, exercising the full encode → channel → decode pipeline "
        "through the public CLI. Results are statistically aggregated per "
        "(profile, SNR, payload) combination. 95 % confidence intervals use "
        "the Wilson score method."
    )
    fig.text(
        0.08,
        0.06,
        textwrap.fill(footnote, width=120),
        ha="left",
        va="bottom",
        fontsize=8,
        color="#666666",
        multialignment="left",
        bbox=dict(
            boxstyle="round,pad=0.5",
            facecolor="#f5f5f5",
            edgecolor="#cccccc",
            linewidth=0.5,
        ),
    )
    return fig


def _page_curves(
    groups: dict,
    profiles: list[str],
    payloads: list[int],
    *,
    metric: str,
    title: str,
    ylabel: str,
    log_scale: bool = False,
    ylim: tuple | None = None,
) -> plt.Figure:
    fig, axes = _grid(len(profiles))
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)

    for ax, profile in zip(axes, profiles):
        ax.set_title(PROFILE_LABELS.get(profile, profile), fontsize=9, pad=4)
        ax.set_xlabel("SNR (dB)", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, which="both" if log_scale else "major", alpha=0.35, linewidth=0.6)
        if log_scale:
            ax.set_yscale("log")
        if ylim:
            ax.set_ylim(*ylim)

        has_data = False
        for pidx, payload in enumerate(payloads):
            snr_map = groups.get(profile, {}).get(payload, {})
            xs, ys = [], []
            for snr in sorted(snr_map):
                st = _stats(snr_map[snr])
                val = st.get(metric, math.nan)
                if _ok(val) and (not log_scale or val > 0):
                    xs.append(snr)
                    ys.append(val)
            if xs:
                has_data = True
                ax.plot(
                    xs,
                    ys,
                    color=_COLORS[pidx % len(_COLORS)],
                    marker=_MARKERS[pidx % len(_MARKERS)],
                    linewidth=1.4,
                    markersize=4,
                    label=f"{payload} bits",
                )

        if has_data:
            ax.legend(fontsize=7, loc="best", framealpha=0.7)
        else:
            ax.text(
                0.5,
                0.5,
                "No data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="gray",
                fontsize=9,
            )

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


def _page_timing(groups: dict, profiles: list[str], payloads: list[int]) -> plt.Figure:
    fig, axes = _grid(len(profiles))
    fig.suptitle(
        "Synchronisation Rate and Timing Error vs. SNR",
        fontsize=13,
        fontweight="bold",
        y=0.99,
    )

    for ax, profile in zip(axes, profiles):
        ax.set_title(PROFILE_LABELS.get(profile, profile), fontsize=9, pad=4)
        ax.set_xlabel("SNR (dB)", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Sync rate", fontsize=8, color="#1f77b4")
        ax.tick_params(labelsize=7, axis="y", labelcolor="#1f77b4")
        ax.grid(True, alpha=0.35, linewidth=0.6)

        ax2 = ax.twinx()
        ax2.set_ylabel("Timing error p95 (samples)", fontsize=7, color="#ff7f0e")
        ax2.tick_params(axis="y", labelcolor="#ff7f0e", labelsize=7)

        for pidx, payload in enumerate(payloads):
            snr_map = groups.get(profile, {}).get(payload, {})
            sxs, sys_, txs, tys = [], [], [], []
            for snr in sorted(snr_map):
                st = _stats(snr_map[snr])
                sr = st.get("sync_rate", math.nan)
                te = st.get("timing_p95", math.nan)
                if _ok(sr):
                    sxs.append(snr)
                    sys_.append(sr)
                if _ok(te):
                    txs.append(snr)
                    tys.append(te)
            if sxs:
                ax.plot(
                    sxs,
                    sys_,
                    color=_COLORS[pidx % len(_COLORS)],
                    marker=_MARKERS[pidx % len(_MARKERS)],
                    linewidth=1.4,
                    markersize=3,
                    label=f"{payload}b",
                )
            if txs:
                ax2.plot(
                    txs,
                    tys,
                    color=_COLORS[pidx % len(_COLORS)],
                    linestyle="--",
                    linewidth=1.0,
                    markersize=2,
                )

        ax.legend(fontsize=7, loc="lower right", framealpha=0.7)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


def _page_summary(records: list[dict[str, Any]], profiles: list[str]) -> plt.Figure:
    fig = plt.figure(figsize=A4)
    fig.suptitle(
        "Campaign Summary — Aggregate Statistics per Profile",
        fontsize=13,
        fontweight="bold",
        y=0.97,
    )

    col_labels = [
        "Profile",
        "Frames",
        "CRC\nrate",
        "FER",
        "BER\n(agg)",
        "Sync\nrate",
        "RT p50\n(ms)",
        "RT p95\n(ms)",
    ]
    rows_data = []
    for profile in profiles:
        subset = [r for r in records if r["profile"] == profile]
        if not subset:
            continue
        st = _stats(subset)
        rows_data.append(
            [
                PROFILE_LABELS.get(profile, profile),
                str(st["n"]),
                f"{st['crc_rate']:.4f}",
                f"{st['fer']:.4f}",
                f"{st['ber_agg']:.2e}" if _ok(st.get("ber_agg", math.nan)) else "—",
                f"{st['sync_rate']:.4f}" if _ok(st.get("sync_rate", math.nan)) else "—",
                f"{st['runtime_p50']:.1f}"
                if _ok(st.get("runtime_p50", math.nan))
                else "—",
                f"{st['runtime_p95']:.1f}"
                if _ok(st.get("runtime_p95", math.nan))
                else "—",
            ]
        )

    ax = fig.add_axes([0.04, 0.20, 0.92, 0.70])
    ax.axis("off")

    if rows_data:
        tbl = ax.table(
            cellText=rows_data,
            colLabels=col_labels,
            loc="center",
            cellLoc="center",
            colWidths=[0.26, 0.09, 0.10, 0.08, 0.12, 0.10, 0.12, 0.13],
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        tbl.scale(1.0, 1.85)
        for (row, col), cell in tbl.get_celld().items():
            cell.set_edgecolor("#bdc3c7")
            cell.get_text().set_wrap(True)
            if col == 0:
                cell.get_text().set_ha("left")
            if row == 0:
                cell.set_facecolor("#2c3e50")
                cell.set_text_props(color="white", fontweight="bold")
            elif row % 2 == 0:
                cell.set_facecolor("#ecf0f1")
            else:
                cell.set_facecolor("white")
    else:
        ax.text(0.5, 0.5, "No data available.", ha="center", va="center", fontsize=12)

    fig.text(
        0.04,
        0.07,
        "Statistics aggregated across all SNR points and payload sizes per profile.",
        fontsize=8,
        color="#666666",
    )
    return fig


# ── Report assembly ───────────────────────────────────────────────────────────


def generate_report(
    records: list[dict[str, Any]],
    output: Path,
    manifest: dict[str, Any],
    version: str,
    title: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    groups = _group(records)

    found = set(groups)
    profiles = [p for p in PROFILE_ORDER if p in found] + sorted(
        found - set(PROFILE_ORDER)
    )
    payloads = sorted({r["payload_bits"] for r in records if r["payload_bits"] > 0})

    if not profiles or not payloads:
        print(
            "[report] WARNING: no groupable records; PDF will be sparse.",
            file=sys.stderr,
        )

    with PdfPages(output) as pdf:
        info = pdf.infodict()
        info["Title"] = title
        info["Author"] = "CDSS validation pipeline"
        info["Subject"] = "Monte Carlo BER/FER campaign"
        info["Creator"] = "generate_report.py"
        info["CreationDate"] = datetime.now(timezone.utc)

        def _save(fig: plt.Figure) -> None:
            # Keep a fixed A4 media box for every page; "tight" can change page size
            # based on text extents and break layout consistency on the cover page.
            pdf.savefig(fig)
            plt.close(fig)

        _save(_page_title(manifest, version, title))

        if profiles and payloads:
            _save(
                _page_curves(
                    groups,
                    profiles,
                    payloads,
                    metric="ber_agg",
                    title="Bit Error Rate vs. SNR",
                    ylabel="BER (aggregate)",
                    log_scale=True,
                    ylim=(1e-5, 1.5),
                )
            )
            _save(
                _page_curves(
                    groups,
                    profiles,
                    payloads,
                    metric="fer",
                    title="Frame Error Rate vs. SNR",
                    ylabel="FER",
                    ylim=(-0.05, 1.05),
                )
            )
            _save(
                _page_curves(
                    groups,
                    profiles,
                    payloads,
                    metric="crc_rate",
                    title="CRC Success Rate vs. SNR",
                    ylabel="CRC success rate",
                    ylim=(-0.05, 1.05),
                )
            )
            _save(_page_timing(groups, profiles, payloads))

        _save(_page_summary(records, profiles))

    size_kb = output.stat().st_size // 1024
    print(
        f"[report] Written {output}  ({size_kb} KB, {len(records)} records, "
        f"{len(profiles)} profiles, {len(payloads)} payload sizes)"
    )


# ── CLI ───────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", required=True, type=Path, metavar="JSONL")
    parser.add_argument("--output", required=True, type=Path, metavar="PDF")
    parser.add_argument("--manifest", default=None, type=Path, metavar="JSON")
    parser.add_argument("--version", default="", metavar="VER")
    parser.add_argument(
        "--title", default="CDSS Modem Validation Report", metavar="TEXT"
    )
    args = parser.parse_args()

    if not args.input.is_file():
        sys.exit(f"[report] ERROR: input not found: {args.input}")

    manifest_path = args.manifest or (args.input.parent / "manifest.json")
    manifest: dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            print(f"[report] WARNING: could not parse manifest: {exc}", file=sys.stderr)

    records = load_records(args.input)
    print(f"[report] Loaded {len(records)} records.", file=sys.stderr)

    generate_report(records, args.output, manifest, args.version, args.title)


if __name__ == "__main__":
    main()
