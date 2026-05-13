#!/usr/bin/env python3
"""Analyze and plot CDSS Monte Carlo JSONL campaign results."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Iterable

RESEARCH_ROOT = Path(__file__).resolve().parents[1]

SUMMARY_FIELDS = [
    "frames",
    "command_failures",
    "decode_attempt_rate",
    "crc_successes",
    "crc_rate",
    "crc_ci95_lo",
    "crc_ci95_hi",
    "fer",
    "total_bits",
    "total_bit_errors",
    "ber_aggregate",
    "ber_ci95_lo",
    "ber_ci95_hi",
    "ber_frame_mean",
    "ber_frame_median",
    "ber_frame_p95",
    "sync_frames",
    "sync_successes",
    "sync_rate",
    "timing_abs_mean",
    "timing_abs_p95",
    "acquisition_score_mean",
    "acquisition_score_p05",
    "acquisition_score_p50",
    "acquisition_score_p95",
    "runtime_ms_mean",
    "runtime_ms_p50",
    "runtime_ms_p95",
    "iterations_mean",
    "iterations_p95",
]


def to_float(value: Any, default: float = math.nan) -> float:
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int(value: Any, default: int = 0) -> int:
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default


def to_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "on"}
    return bool(value)


def finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def mean(values: Iterable[float]) -> float:
    vals = [float(v) for v in values if finite(v)]
    return sum(vals) / len(vals) if vals else math.nan


def percentile(values: Iterable[float], p: float) -> float:
    vals = sorted(float(v) for v in values if finite(v))
    if not vals:
        return math.nan
    if p <= 0.0:
        return vals[0]
    if p >= 100.0:
        return vals[-1]
    rank = (len(vals) - 1) * p / 100.0
    lo = math.floor(rank)
    hi = math.ceil(rank)
    if lo == hi:
        return vals[int(lo)]
    frac = rank - lo
    return vals[int(lo)] * (1.0 - frac) + vals[int(hi)] * frac


def wilson_interval(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return math.nan, math.nan
    phat = successes / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2.0 * n)) / denom
    margin = (z / denom) * math.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n))
    return max(0.0, center - margin), min(1.0, center + margin)


def fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "nan"
    if isinstance(value, float) and not math.isfinite(value):
        return "nan"
    if isinstance(value, (int, float)):
        return f"{float(value):.{digits}f}"
    return str(value)


def slug(text: Any) -> str:
    chars = []
    for char in str(text).lower():
        if char.isalnum():
            chars.append(char)
        elif char in {"-", "_"}:
            chars.append(char)
        else:
            chars.append("_")
    return "".join(chars).strip("_") or "unknown"


def normalize_record(raw: dict[str, Any], row_number: int) -> dict[str, Any]:
    experiment = raw.get("experiment", {}) or {}
    signal = raw.get("signal", {}) or {}
    channel = raw.get("channel", {}) or {}
    acquisition = raw.get("acquisition", {}) or {}
    decoder = raw.get("decoder", {}) or {}
    payload = raw.get("payload", {}) or {}
    performance = raw.get("performance", {}) or {}
    process = raw.get("process", {}) or {}

    profile = str(
        experiment.get("profile")
        or channel.get("profile")
        or raw.get("profile")
        or "legacy"
    )
    payload_bits = to_int(signal.get("payload_bits", raw.get("payload_bits", payload.get("total_bits", -1))), -1)
    snr_db = to_float(channel.get("snr_db", raw.get("snr_db")))
    exit_code = to_int(experiment.get("exit_code", raw.get("exit_code", 0)), 0)
    parse_error = str(process.get("parse_error", raw.get("parse_error", "")) or "")

    bit_errors = to_float(payload.get("bit_errors", raw.get("bit_errors")))
    total_bits = to_float(payload.get("total_bits", payload_bits if payload_bits > 0 else math.nan))
    ber = to_float(payload.get("ber", payload.get("ber_frame", raw.get("ber"))))
    if not finite(ber) and finite(bit_errors) and finite(total_bits) and total_bits > 0:
        ber = bit_errors / total_bits
    if not finite(bit_errors) and finite(ber) and finite(total_bits):
        bit_errors = round(ber * total_bits)

    sync_search = to_bool(acquisition.get("sync_search"), False)
    frame_detected = to_bool(acquisition.get("frame_detected"), False)
    sync_ok = to_bool(acquisition.get("sync_ok"), frame_detected if sync_search else True)

    iterations = to_float(
        decoder.get(
            "refinement_iters_used",
            decoder.get("iterations_used", decoder.get("ldpc_iterations")),
        )
    )

    command = experiment.get("command", raw.get("command", ""))
    if isinstance(command, list):
        command_text = " ".join(str(part) for part in command)
    else:
        command_text = str(command)

    return {
        "row": row_number,
        "profile": profile,
        "trial": to_int(experiment.get("trial", raw.get("trial", -1)), -1),
        "seed": to_int(experiment.get("seed", raw.get("seed", -1)), -1),
        "payload_bits": payload_bits,
        "snr_db": snr_db,
        "exit_code": exit_code,
        "command": command_text,
        "parse_error": parse_error,
        "stderr": str(process.get("stderr", raw.get("stderr", "")) or ""),
        "runtime_ms": to_float(experiment.get("runtime_ms", performance.get("runtime_ms"))),
        "sync_search": sync_search,
        "frame_detected": frame_detected,
        "sync_ok": sync_ok,
        "acquisition_score": to_float(acquisition.get("acquisition_score")),
        "timing_error_samples": to_float(acquisition.get("timing_error_samples")),
        "decode_attempted": to_bool(decoder.get("decode_attempted"), False),
        "crc_ok": to_bool(decoder.get("crc_ok"), False),
        "iterations_used": iterations,
        "bit_errors": bit_errors,
        "total_bits": total_bits,
        "ber": ber,
    }


def load_records(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    parse_failures: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as stream:
        for row_number, line in enumerate(stream, 1):
            text = line.strip()
            if not text:
                continue
            try:
                raw = json.loads(text)
            except json.JSONDecodeError as exc:
                parse_failures.append({
                    "source": "input_jsonl",
                    "row": row_number,
                    "profile": "",
                    "payload_bits": "",
                    "snr_db": "",
                    "trial": "",
                    "seed": "",
                    "exit_code": "",
                    "parse_error": str(exc),
                    "stderr": "",
                })
                continue
            records.append(normalize_record(raw, row_number))
    return records, parse_failures


def summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    if n == 0:
        return {field: math.nan for field in SUMMARY_FIELDS}

    command_failures = sum(1 for r in records if r["exit_code"] != 0 or r["parse_error"])
    decode_attempts = sum(1 for r in records if r["decode_attempted"])
    crc_successes = sum(1 for r in records if r["crc_ok"])
    crc_lo, crc_hi = wilson_interval(crc_successes, n)

    bit_error_values = [r["bit_errors"] for r in records if finite(r["bit_errors"])]
    total_bit_values = [r["total_bits"] for r in records if finite(r["total_bits"])]
    total_bit_errors = int(sum(bit_error_values)) if bit_error_values else 0
    total_bits = int(sum(total_bit_values)) if total_bit_values else 0
    ber_lo, ber_hi = wilson_interval(total_bit_errors, total_bits)

    ber_frames = [r["ber"] for r in records if finite(r["ber"])]
    sync_records = [r for r in records if r["sync_search"]]
    sync_successes = sum(1 for r in sync_records if r["sync_ok"])
    timing_abs = [abs(r["timing_error_samples"]) for r in sync_records if finite(r["timing_error_samples"])]
    acquisition_scores = [r["acquisition_score"] for r in sync_records if finite(r["acquisition_score"])]
    runtimes = [r["runtime_ms"] for r in records if finite(r["runtime_ms"])]
    iterations = [r["iterations_used"] for r in records if finite(r["iterations_used"])]

    return {
        "frames": n,
        "command_failures": command_failures,
        "decode_attempt_rate": decode_attempts / n,
        "crc_successes": crc_successes,
        "crc_rate": crc_successes / n,
        "crc_ci95_lo": crc_lo,
        "crc_ci95_hi": crc_hi,
        "fer": 1.0 - crc_successes / n,
        "total_bits": total_bits,
        "total_bit_errors": total_bit_errors,
        "ber_aggregate": total_bit_errors / total_bits if total_bits > 0 else math.nan,
        "ber_ci95_lo": ber_lo,
        "ber_ci95_hi": ber_hi,
        "ber_frame_mean": mean(ber_frames),
        "ber_frame_median": median(ber_frames) if ber_frames else math.nan,
        "ber_frame_p95": percentile(ber_frames, 95.0),
        "sync_frames": len(sync_records),
        "sync_successes": sync_successes,
        "sync_rate": sync_successes / len(sync_records) if sync_records else math.nan,
        "timing_abs_mean": mean(timing_abs),
        "timing_abs_p95": percentile(timing_abs, 95.0),
        "acquisition_score_mean": mean(acquisition_scores),
        "acquisition_score_p05": percentile(acquisition_scores, 5.0),
        "acquisition_score_p50": percentile(acquisition_scores, 50.0),
        "acquisition_score_p95": percentile(acquisition_scores, 95.0),
        "runtime_ms_mean": mean(runtimes),
        "runtime_ms_p50": percentile(runtimes, 50.0),
        "runtime_ms_p95": percentile(runtimes, 95.0),
        "iterations_mean": mean(iterations),
        "iterations_p95": percentile(iterations, 95.0),
    }


def group_records(records: list[dict[str, Any]], keys: tuple[str, ...]) -> dict[tuple[Any, ...], list[dict[str, Any]]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        values: list[Any] = []
        skip = False
        for key in keys:
            value = record[key]
            if key in {"payload_bits", "snr_db"} and not finite(value):
                skip = True
                break
            values.append(value)
        if not skip:
            grouped[tuple(values)].append(record)
    return dict(grouped)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def snr_threshold(points: list[tuple[float, float]], target: float) -> float:
    valid_points = sorted((float(s), float(r)) for s, r in points if finite(s) and finite(r))
    if not valid_points:
        return math.nan
    if all(rate >= target for _, rate in valid_points):
        return valid_points[0][0]
    for idx in range(1, len(valid_points)):
        s0, r0 = valid_points[idx - 1]
        s1, r1 = valid_points[idx]
        if (r0 <= target <= r1) or (r1 <= target <= r0):
            if r0 == r1:
                return s1
            return s0 + (target - r0) * (s1 - s0) / (r1 - r0)
    return math.nan


def load_manifest(input_path: Path, explicit_path: Path | None) -> dict[str, Any]:
    candidates = []
    if explicit_path is not None:
        candidates.append(explicit_path)
    candidates.append(input_path.parent / "manifest.json")
    for candidate in candidates:
        if candidate.exists():
            try:
                return json.loads(candidate.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                return {}
    return {}


def build_summary_rows(records: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key_values, group in sorted(group_records(records, keys).items()):
        row = {key: value for key, value in zip(keys, key_values)}
        row.update(summarize(group))
        rows.append(row)
    return rows


def build_threshold_rows(rows_by_profile_payload_snr: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int], list[tuple[float, float]]] = defaultdict(list)
    for row in rows_by_profile_payload_snr:
        grouped[(str(row["profile"]), int(row["payload_bits"]))].append(
            (float(row["snr_db"]), float(row["crc_rate"]))
        )

    rows: list[dict[str, Any]] = []
    for (profile, payload_bits), points in sorted(grouped.items()):
        rows.append({
            "profile": profile,
            "payload_bits": payload_bits,
            "snr_at_50": snr_threshold(points, 0.50),
            "snr_at_90": snr_threshold(points, 0.90),
            "snr_at_95": snr_threshold(points, 0.95),
            "snr_at_99": snr_threshold(points, 0.99),
        })
    return rows


def build_failure_rows(records: list[dict[str, Any]], parse_failures: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = list(parse_failures)
    for record in records:
        if record["exit_code"] == 0 and not record["parse_error"]:
            continue
        rows.append({
            "source": "trial",
            "row": record["row"],
            "profile": record["profile"],
            "payload_bits": record["payload_bits"],
            "snr_db": record["snr_db"],
            "trial": record["trial"],
            "seed": record["seed"],
            "exit_code": record["exit_code"],
            "parse_error": record["parse_error"],
            "stderr": record["stderr"],
        })
    return rows


def write_report(
    path: Path,
    manifest: dict[str, Any],
    input_path: Path,
    records: list[dict[str, Any]],
    parse_failures: list[dict[str, Any]],
    by_profile_rows: list[dict[str, Any]],
    threshold_rows: list[dict[str, Any]],
    figure_files: list[str],
) -> None:
    global_summary = summarize(records)
    profiles = sorted({str(r["profile"]) for r in records})
    payloads = sorted({int(r["payload_bits"]) for r in records if r["payload_bits"] >= 0})
    snrs = sorted({float(r["snr_db"]) for r in records if finite(r["snr_db"])})

    lines = [
        "# CDSS Monte Carlo Analysis",
        "",
        "## Dataset",
        "",
        f"- Input JSONL: `{input_path}`",
        f"- Manifest created UTC: `{manifest.get('created_utc', 'unknown')}`",
        f"- Binary: `{manifest.get('binary', 'unknown')}`",
        f"- Binary SHA-256: `{manifest.get('binary_sha256', 'unknown')}`",
        f"- Frames analyzed: **{global_summary['frames']}**",
        f"- Input JSON parse failures: **{len(parse_failures)}**",
        f"- Profiles: **{profiles}**",
        f"- Payload sizes (bits): **{payloads}**",
        f"- SNR points (dB): **{snrs}**",
        "",
        "## Global Performance",
        "",
        f"- CRC success: **{fmt(global_summary['crc_rate'] * 100.0, 2)}%** "
        f"(95% CI {fmt(global_summary['crc_ci95_lo'] * 100.0, 2)}% to {fmt(global_summary['crc_ci95_hi'] * 100.0, 2)}%)",
        f"- FER: **{fmt(global_summary['fer'] * 100.0, 2)}%**",
        f"- Aggregate BER: **{fmt(global_summary['ber_aggregate'], 8)}** "
        f"(95% CI {fmt(global_summary['ber_ci95_lo'], 8)} to {fmt(global_summary['ber_ci95_hi'], 8)})",
        f"- Command failures: **{global_summary['command_failures']}**",
        f"- Runtime mean / P95: **{fmt(global_summary['runtime_ms_mean'], 2)} ms / {fmt(global_summary['runtime_ms_p95'], 2)} ms**",
        "",
        "## By Profile",
        "",
        "| Profile | Frames | CRC % | FER % | BER aggregate | Sync % | Runtime mean ms |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in by_profile_rows:
        sync_text = "n/a" if not finite(row["sync_rate"]) else fmt(row["sync_rate"] * 100.0, 2)
        lines.append(
            f"| {row['profile']} | {row['frames']} | {fmt(row['crc_rate'] * 100.0, 2)} | "
            f"{fmt(row['fer'] * 100.0, 2)} | {fmt(row['ber_aggregate'], 8)} | "
            f"{sync_text} | {fmt(row['runtime_ms_mean'], 2)} |"
        )

    lines.extend([
        "",
        "## CRC SNR Thresholds",
        "",
        "| Profile | Payload bits | 50% | 90% | 95% | 99% |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in threshold_rows:
        lines.append(
            f"| {row['profile']} | {row['payload_bits']} | {fmt(row['snr_at_50'], 2)} | "
            f"{fmt(row['snr_at_90'], 2)} | {fmt(row['snr_at_95'], 2)} | {fmt(row['snr_at_99'], 2)} |"
        )

    lines.extend([
        "",
        "## Figures",
        "",
    ])
    if figure_files:
        for figure in figure_files:
            lines.append(f"- `{figure}`")
    else:
        lines.append("- No figures were generated because matplotlib was unavailable or the dataset was empty.")

    lines.extend([
        "",
        "## Interpretation Boundaries",
        "",
        "- These plots and tables summarize the measured `cdss simulate` campaign only.",
        "- The CSV tables are the authoritative numeric artifacts; PNG files are derived visualizations.",
        "- Profiles using synchronized acquisition exercise preamble search and phase alignment before decoding.",
        "- SNR thresholds are linear interpolations over measured SNR points; single-point smoke campaigns are not threshold measurements.",
        "- Scientific claims should be restricted to the measured payload, SNR, profile, and seed grid unless additional campaigns are run.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def import_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except Exception as exc:  # noqa: BLE001 - reported to caller as optional plotting failure.
        print(f"Warning: matplotlib unavailable, skipping PNG figures: {exc}", file=sys.stderr)
        return None


def plot_metric_by_group(
    plt: Any,
    rows: list[dict[str, Any]],
    out_path: Path,
    group_key: str,
    metric_key: str,
    title: str,
    ylabel: str,
    scale: float = 1.0,
    semilogy: bool = False,
    y_floor: float | None = None,
    ylim: tuple[float, float] | None = None,
) -> bool:
    grouped: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if finite(row.get("snr_db")) and finite(row.get(metric_key)):
            grouped[row[group_key]].append(row)
    if not grouped:
        return False

    plt.figure(figsize=(10, 6))
    for group, group_rows in sorted(grouped.items(), key=lambda item: str(item[0])):
        points = sorted(group_rows, key=lambda row: float(row["snr_db"]))
        xs = [float(row["snr_db"]) for row in points]
        ys = [float(row[metric_key]) * scale for row in points]
        if y_floor is not None:
            ys = [max(y, y_floor) for y in ys]
        if semilogy:
            plt.semilogy(xs, ys, marker="o", linewidth=1.5, markersize=4, label=str(group))
        else:
            plt.plot(xs, ys, marker="o", linewidth=1.5, markersize=4, label=str(group))
    plt.title(title)
    plt.xlabel("SNR (dB)")
    plt.ylabel(ylabel)
    if ylim is not None and not semilogy:
        plt.ylim(*ylim)
    plt.grid(True, which="both" if semilogy else "major", linestyle="--", alpha=0.4)
    plt.legend(title=group_key.replace("_", " ").title())
    plt.tight_layout()
    plt.savefig(out_path, dpi=250)
    plt.close()
    return True


def plot_runtime(plt: Any, rows: list[dict[str, Any]], out_path: Path) -> bool:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if finite(row.get("snr_db")):
            grouped[str(row["profile"])].append(row)
    if not grouped:
        return False

    plt.figure(figsize=(10, 6))
    for profile, group_rows in sorted(grouped.items()):
        points = sorted(group_rows, key=lambda row: float(row["snr_db"]))
        xs = [float(row["snr_db"]) for row in points]
        p50 = [float(row["runtime_ms_p50"]) for row in points]
        p95 = [float(row["runtime_ms_p95"]) for row in points]
        plt.plot(xs, p50, marker="o", linewidth=1.5, markersize=4, label=f"{profile} p50")
        plt.plot(xs, p95, marker="s", linewidth=1.1, markersize=3, linestyle="--", label=f"{profile} p95")
    plt.title("Runtime vs SNR by Profile")
    plt.xlabel("SNR (dB)")
    plt.ylabel("Runtime (ms)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(fontsize="small", ncol=2)
    plt.tight_layout()
    plt.savefig(out_path, dpi=250)
    plt.close()
    return True


def plot_thresholds(plt: Any, rows: list[dict[str, Any]], out_path: Path) -> bool:
    valid_rows = [row for row in rows if finite(row.get("snr_at_95"))]
    if not valid_rows:
        return False
    labels = [f"{row['profile']}\n{row['payload_bits']}b" for row in valid_rows]
    values = [float(row["snr_at_95"]) for row in valid_rows]
    width = max(10, min(22, 0.45 * len(labels) + 4))
    plt.figure(figsize=(width, 6))
    plt.bar(range(len(labels)), values)
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.ylabel("SNR at 95% CRC success (dB)")
    plt.title("Estimated 95% CRC Success Threshold")
    plt.grid(True, axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=250)
    plt.close()
    return True


def plot_heatmaps(plt: Any, rows: list[dict[str, Any]], out_dir: Path) -> list[str]:
    figures: list[str] = []
    profiles = sorted({str(row["profile"]) for row in rows})
    for profile in profiles:
        profile_rows = [row for row in rows if str(row["profile"]) == profile]
        payloads = sorted({int(row["payload_bits"]) for row in profile_rows})
        snrs = sorted({float(row["snr_db"]) for row in profile_rows})
        if not payloads or not snrs:
            continue
        index = {
            (int(row["payload_bits"]), float(row["snr_db"])): row
            for row in profile_rows
        }
        heatmap_specs = [
            ("crc_rate", "CRC Success (%)", 100.0, None, f"heatmap_crc_success_{slug(profile)}.png"),
            ("ber_aggregate", "log10 Aggregate BER", 1.0, "log10", f"heatmap_log10_ber_{slug(profile)}.png"),
        ]
        for metric_key, color_label, scale, transform, filename in heatmap_specs:
            matrix: list[list[float]] = []
            for payload_bits in payloads:
                row_values = []
                for snr_db in snrs:
                    value = index.get((payload_bits, snr_db), {}).get(metric_key, math.nan)
                    if transform == "log10":
                        value = math.log10(max(float(value), 1e-12)) if finite(value) else math.nan
                    elif finite(value):
                        value = float(value) * scale
                    row_values.append(value)
                matrix.append(row_values)
            plt.figure(figsize=(10, max(4, 0.6 * len(payloads) + 2)))
            image = plt.imshow(matrix, aspect="auto", origin="lower")
            plt.xticks(range(len(snrs)), [fmt(s, 1) for s in snrs], rotation=45)
            plt.yticks(range(len(payloads)), [str(p) for p in payloads])
            plt.xlabel("SNR (dB)")
            plt.ylabel("Payload bits")
            plt.title(f"{color_label} Heatmap - {profile}")
            plt.colorbar(image, label=color_label)
            plt.tight_layout()
            plt.savefig(out_dir / filename, dpi=250)
            plt.close()
            figures.append(filename)
    return figures


def generate_figures(
    out_dir: Path,
    by_profile_snr_rows: list[dict[str, Any]],
    by_profile_payload_snr_rows: list[dict[str, Any]],
    threshold_rows: list[dict[str, Any]],
) -> list[str]:
    plt = import_matplotlib()
    if plt is None:
        return []

    figures: list[str] = []
    specs = [
        (
            by_profile_snr_rows,
            "crc_success_vs_snr_by_profile.png",
            "profile",
            "crc_rate",
            "CRC Success vs SNR by Profile",
            "CRC Success (%)",
            100.0,
            False,
            None,
            (-2.0, 102.0),
        ),
        (
            by_profile_snr_rows,
            "fer_vs_snr_by_profile.png",
            "profile",
            "fer",
            "Frame Error Rate vs SNR by Profile",
            "FER (%)",
            100.0,
            False,
            None,
            (-2.0, 102.0),
        ),
        (
            by_profile_snr_rows,
            "ber_vs_snr_by_profile.png",
            "profile",
            "ber_aggregate",
            "Aggregate BER vs SNR by Profile",
            "Aggregate BER",
            1.0,
            True,
            1e-12,
            None,
        ),
        (
            by_profile_snr_rows,
            "sync_success_vs_snr_by_profile.png",
            "profile",
            "sync_rate",
            "Synchronization Success vs SNR by Profile",
            "Sync Success (%)",
            100.0,
            False,
            None,
            (-2.0, 102.0),
        ),
        (
            by_profile_snr_rows,
            "timing_error_p95_vs_snr_by_profile.png",
            "profile",
            "timing_abs_p95",
            "P95 Absolute Timing Error vs SNR by Profile",
            "P95 absolute timing error (samples)",
            1.0,
            False,
            None,
            None,
        ),
        (
            by_profile_snr_rows,
            "acquisition_score_vs_snr_by_profile.png",
            "profile",
            "acquisition_score_mean",
            "Mean Acquisition Score vs SNR by Profile",
            "Mean acquisition score",
            1.0,
            False,
            None,
            None,
        ),
    ]
    for rows, filename, group_key, metric_key, title, ylabel, scale, semilogy, y_floor, ylim in specs:
        if plot_metric_by_group(
            plt,
            rows,
            out_dir / filename,
            group_key,
            metric_key,
            title,
            ylabel,
            scale=scale,
            semilogy=semilogy,
            y_floor=y_floor,
            ylim=ylim,
        ):
            figures.append(filename)

    if plot_runtime(plt, by_profile_snr_rows, out_dir / "runtime_ms_vs_snr_by_profile.png"):
        figures.append("runtime_ms_vs_snr_by_profile.png")

    profiles = sorted({str(row["profile"]) for row in by_profile_payload_snr_rows})
    for profile in profiles:
        rows = [row for row in by_profile_payload_snr_rows if str(row["profile"]) == profile]
        crc_file = f"crc_success_vs_snr_{slug(profile)}_by_payload.png"
        if plot_metric_by_group(
            plt,
            rows,
            out_dir / crc_file,
            "payload_bits",
            "crc_rate",
            f"CRC Success vs SNR - {profile}",
            "CRC Success (%)",
            scale=100.0,
            ylim=(-2.0, 102.0),
        ):
            figures.append(crc_file)
        ber_file = f"ber_vs_snr_{slug(profile)}_by_payload.png"
        if plot_metric_by_group(
            plt,
            rows,
            out_dir / ber_file,
            "payload_bits",
            "ber_aggregate",
            f"Aggregate BER vs SNR - {profile}",
            "Aggregate BER",
            semilogy=True,
            y_floor=1e-12,
        ):
            figures.append(ber_file)

    heatmap_figures = plot_heatmaps(plt, by_profile_payload_snr_rows, out_dir)
    figures.extend(heatmap_figures)

    threshold_file = "snr_threshold_crc95_by_profile_payload.png"
    if plot_thresholds(plt, threshold_rows, out_dir / threshold_file):
        figures.append(threshold_file)

    return figures


def default_input_path() -> Path:
    last_run = RESEARCH_ROOT / ".last_run"
    if last_run.exists():
        run_dir_text = last_run.read_text(encoding="utf-8").strip()
        if run_dir_text:
            run_dir = Path(run_dir_text)
            if not run_dir.is_absolute():
                run_dir = RESEARCH_ROOT / run_dir
            return run_dir / "results.jsonl"
    return RESEARCH_ROOT / "mc_results" / "cdss_data.jsonl"


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze and plot CDSS Monte Carlo JSONL results.")
    parser.add_argument("--input", help="Input results.jsonl file. Defaults to research/.last_run when available.")
    parser.add_argument("--out-dir", help="Output directory. Defaults to <input directory>/analysis.")
    parser.add_argument("--manifest", help="Optional manifest.json path. Defaults to <input directory>/manifest.json.")
    args = parser.parse_args()

    input_path = Path(args.input).resolve() if args.input else default_input_path().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input JSONL file not found: {input_path}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else (input_path.parent / "analysis").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = Path(args.manifest).resolve() if args.manifest else None
    manifest = load_manifest(input_path, manifest_path)
    records, parse_failures = load_records(input_path)
    if not records:
        raise RuntimeError("No Monte Carlo records found in input JSONL file")

    global_row = {"scope": "all", "input_parse_failures": len(parse_failures)}
    global_row.update(summarize(records))
    by_profile_payload_snr_rows = build_summary_rows(records, ("profile", "payload_bits", "snr_db"))
    by_profile_payload_rows = build_summary_rows(records, ("profile", "payload_bits"))
    by_profile_snr_rows = build_summary_rows(records, ("profile", "snr_db"))
    by_profile_rows = build_summary_rows(records, ("profile",))
    threshold_rows = build_threshold_rows(by_profile_payload_snr_rows)
    failure_rows = build_failure_rows(records, parse_failures)

    write_csv(out_dir / "global_summary.csv", [global_row], ["scope", "input_parse_failures", *SUMMARY_FIELDS])
    write_csv(
        out_dir / "summary_by_profile_payload_snr.csv",
        by_profile_payload_snr_rows,
        ["profile", "payload_bits", "snr_db", *SUMMARY_FIELDS],
    )
    write_csv(
        out_dir / "summary_by_profile_payload.csv",
        by_profile_payload_rows,
        ["profile", "payload_bits", *SUMMARY_FIELDS],
    )
    write_csv(
        out_dir / "summary_by_profile_snr.csv",
        by_profile_snr_rows,
        ["profile", "snr_db", *SUMMARY_FIELDS],
    )
    write_csv(
        out_dir / "summary_by_profile.csv",
        by_profile_rows,
        ["profile", *SUMMARY_FIELDS],
    )
    write_csv(
        out_dir / "snr_thresholds_crc.csv",
        threshold_rows,
        ["profile", "payload_bits", "snr_at_50", "snr_at_90", "snr_at_95", "snr_at_99"],
    )
    write_csv(
        out_dir / "failures.csv",
        failure_rows,
        ["source", "row", "profile", "payload_bits", "snr_db", "trial", "seed", "exit_code", "parse_error", "stderr"],
    )

    figure_files = generate_figures(out_dir, by_profile_snr_rows, by_profile_payload_snr_rows, threshold_rows)
    write_report(
        out_dir / "analysis_summary.md",
        manifest,
        input_path,
        records,
        parse_failures,
        by_profile_rows,
        threshold_rows,
        figure_files,
    )

    print(f"Analysis complete. Tables and figures written to: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
