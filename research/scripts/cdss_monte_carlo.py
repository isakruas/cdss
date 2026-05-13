#!/usr/bin/env python3
"""Run a reproducible CDSS Monte Carlo campaign through `cdss simulate`.

The campaign intentionally exercises the public command-line modem path instead
of reimplementing the channel in Python. Each trial calls:

    cdss simulate --bits N --snr-db X --seed S [profile options]

Outputs are written as JSONL plus CSV/Markdown summaries suitable for review,
regression tracking, and scientific reports.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_ROOT = REPO_ROOT / "research" / "mc_results"

PROFILE_OPTIONS: dict[str, list[str]] = {
    "aligned_awgn": [],
    "sync_awgn": [
        "--sync-search",
        "--prefix-pad-s",
        "0.020",
        "--suffix-pad-s",
        "0.010",
    ],
    "cfo_drift": [
        "--sync-search",
        "--prefix-pad-s",
        "0.020",
        "--suffix-pad-s",
        "0.010",
        "--cfo-hz",
        "0.020",
        "--drift-hz-s",
        "0.001",
    ],
    "tone": [
        "--sync-search",
        "--prefix-pad-s",
        "0.020",
        "--suffix-pad-s",
        "0.010",
        "--tone-frequency-hz",
        "1525.0",
        "--tone-sir-db",
        "25.0",
    ],
    "clock": [
        "--sync-search",
        "--prefix-pad-s",
        "0.020",
        "--suffix-pad-s",
        "0.010",
        "--clock-ppm",
        "5.0",
    ],
}

PRESETS: dict[str, dict[str, Any]] = {
    "smoke": {
        "payloads": [16],
        "snrs": [45.0],
        "profiles": ["aligned_awgn", "sync_awgn"],
        "trials_per_point": 2,
    },
    "quick": {
        "payloads": [16, 64],
        "snrs": [45.0, 35.0, 25.0],
        "profiles": ["aligned_awgn", "sync_awgn", "cfo_drift", "tone", "clock"],
        "trials_per_point": 10,
    },
    "publication": {
        "payloads": [32, 64, 128, 256],
        "snrs": [45.0, 40.0, 35.0, 30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0,
                 -5.0, -10.0, -15.0, -20.0, -25.0, -30.0, -35.0],
        "profiles": ["aligned_awgn", "sync_awgn", "cfo_drift", "tone", "clock"],
        "trials_per_point": 200,
    },
}


@dataclass(frozen=True)
class Job:
    job_id: int
    profile: str
    payload_bits: int
    snr_db: float
    trial: int
    seed: int


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def parse_int_list(text: str) -> list[int]:
    return [int(item.strip()) for item in text.split(",") if item.strip()]


def parse_float_list(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def parse_str_list(text: str) -> list[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def fmt_float(value: float) -> str:
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.6g}"


def newest_built_binary() -> Path | None:
    candidates = list((REPO_ROOT / "build").glob("**/app/cdss"))
    candidates = [p for p in candidates if p.is_file() and os.access(p, os.X_OK)]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def resolve_binary(given: str | None) -> Path:
    if given:
        path = Path(given)
        if path.exists():
            return path.resolve()
        found = shutil.which(given)
        if found:
            return Path(found).resolve()
        raise FileNotFoundError(f"CDSS binary not found: {given}")

    built = newest_built_binary()
    if built:
        return built.resolve()

    found = shutil.which("cdss")
    if found:
        return Path(found).resolve()

    raise FileNotFoundError(
        "No cdss executable found. Run `make build` or pass --binary /path/to/cdss."
    )


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def first_json_object(stdout: str) -> dict[str, Any]:
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)
    raise ValueError("No JSON object found in cdss simulate output")


def run_trial(binary: Path, job: Job, timeout_s: float, omp_threads: int) -> dict[str, Any]:
    args = [
        str(binary),
        "simulate",
        "--bits",
        str(job.payload_bits),
        "--snr-db",
        fmt_float(job.snr_db),
        "--seed",
        str(job.seed),
        *PROFILE_OPTIONS[job.profile],
    ]
    t0 = time.perf_counter()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(max(1, omp_threads))
    proc = subprocess.run(args, capture_output=True, text=True, timeout=timeout_s, env=env)
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    parsed: dict[str, Any] = {}
    parse_error = ""
    if proc.returncode == 0:
        try:
            parsed = first_json_object(proc.stdout)
        except Exception as exc:  # noqa: BLE001 - stored in result for diagnosis.
            parse_error = str(exc)

    crc_ok = bool(parsed.get("crc_ok", False)) if parsed else False
    bit_errors = int(parsed.get("bit_errors", job.payload_bits)) if parsed else job.payload_bits
    ber = float(parsed.get("ber", bit_errors / max(job.payload_bits, 1))) if parsed else 1.0
    payload_bit_count = int(parsed.get("payload_bit_count", 0)) if parsed else 0
    sync_search = bool(parsed.get("sync_search", "--sync-search" in PROFILE_OPTIONS[job.profile]))
    sync_ok = bool(parsed.get("sync_ok", False)) if sync_search else True

    return {
        "schema_version": "cdss-simulate-monte-carlo-v1",
        "experiment": {
            "job_id": job.job_id,
            "profile": job.profile,
            "trial": job.trial,
            "seed": job.seed,
            "timestamp_utc": utc_now(),
            "command": args,
            "exit_code": proc.returncode,
            "runtime_ms": elapsed_ms,
        },
        "signal": {
            "payload_bits": job.payload_bits,
        },
        "channel": {
            "snr_db": job.snr_db,
            "cfo_hz": float(parsed.get("cfo_hz", 0.0)) if parsed else 0.0,
            "drift_hz_s": float(parsed.get("drift_hz_s", 0.0)) if parsed else 0.0,
            "clock_ppm": float(parsed.get("clock_ppm", 0.0)) if parsed else 0.0,
            "prefix_pad_s": float(parsed.get("prefix_pad_s", 0.0)) if parsed else 0.0,
            "suffix_pad_s": float(parsed.get("suffix_pad_s", 0.0)) if parsed else 0.0,
            "tone_count": int(parsed.get("tone_count", 0)) if parsed else 0,
            "tone_frequency_hz": parsed.get("tone_frequency_hz") if parsed else None,
            "tone_sir_db": parsed.get("tone_sir_db") if parsed else None,
            "noise_var": parsed.get("noise_var") if parsed else None,
        },
        "acquisition": {
            "sync_search": sync_search,
            "sync_ok": sync_ok,
            "expected_sample": parsed.get("expected_sample") if parsed else None,
            "detected_sample": parsed.get("detected_sample") if parsed else None,
            "timing_error_samples": parsed.get("timing_error_samples") if parsed else None,
            "acquisition_score": parsed.get("acquisition_score") if parsed else None,
        },
        "decoder": {
            "decode_attempted": proc.returncode == 0 and bool(parsed),
            "crc_ok": crc_ok,
            "refinement_iters_used": int(parsed.get("refinement_iters_used", 0)) if parsed else 0,
            "payload_bit_count": payload_bit_count,
        },
        "payload": {
            "bit_errors": bit_errors,
            "total_bits": job.payload_bits,
            "ber": ber,
        },
        "process": {
            "stdout": proc.stdout.strip(),
            "stderr": proc.stderr.strip(),
            "parse_error": parse_error,
        },
    }


def wilson_interval(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return (math.nan, math.nan)
    phat = successes / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2.0 * n)) / denom
    margin = (z / denom) * math.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n))
    return max(0.0, center - margin), min(1.0, center + margin)


def mean(values: Iterable[float]) -> float:
    vals = list(values)
    return sum(vals) / len(vals) if vals else math.nan


def percentile(values: Iterable[float], p: float) -> float:
    vals = sorted(values)
    if not vals:
        return math.nan
    if p <= 0:
        return vals[0]
    if p >= 100:
        return vals[-1]
    rank = (len(vals) - 1) * p / 100.0
    lo = math.floor(rank)
    hi = math.ceil(rank)
    if lo == hi:
        return vals[int(lo)]
    frac = rank - lo
    return vals[int(lo)] * (1.0 - frac) + vals[int(hi)] * frac


def summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    crc_successes = sum(1 for r in records if r["decoder"]["crc_ok"])
    command_failures = sum(1 for r in records if r["experiment"]["exit_code"] != 0)
    decode_attempts = sum(1 for r in records if r["decoder"]["decode_attempted"])
    sync_records = [r for r in records if r["acquisition"].get("sync_search")]
    sync_successes = sum(1 for r in sync_records if r["acquisition"].get("sync_ok"))
    total_bits = sum(int(r["payload"]["total_bits"]) for r in records)
    total_bit_errors = sum(int(r["payload"]["bit_errors"]) for r in records)
    crc_lo, crc_hi = wilson_interval(crc_successes, n)
    ber_lo, ber_hi = wilson_interval(total_bit_errors, total_bits)
    runtimes = [float(r["experiment"]["runtime_ms"]) for r in records]
    bers = [float(r["payload"]["ber"]) for r in records]
    timing_abs = [
        abs(float(r["acquisition"]["timing_error_samples"]))
        for r in sync_records
        if r["acquisition"].get("timing_error_samples") is not None
    ]
    acq_scores = [
        float(r["acquisition"]["acquisition_score"])
        for r in sync_records
        if r["acquisition"].get("acquisition_score") is not None
    ]
    return {
        "frames": n,
        "command_failures": command_failures,
        "decode_attempt_rate": decode_attempts / n if n else math.nan,
        "crc_successes": crc_successes,
        "crc_rate": crc_successes / n if n else math.nan,
        "crc_ci95_lo": crc_lo,
        "crc_ci95_hi": crc_hi,
        "fer": 1.0 - crc_successes / n if n else math.nan,
        "total_bits": total_bits,
        "total_bit_errors": total_bit_errors,
        "ber_aggregate": total_bit_errors / total_bits if total_bits else math.nan,
        "ber_ci95_lo": ber_lo,
        "ber_ci95_hi": ber_hi,
        "ber_frame_mean": mean(bers),
        "ber_frame_median": median(bers) if bers else math.nan,
        "ber_frame_p95": percentile(bers, 95),
        "sync_frames": len(sync_records),
        "sync_successes": sync_successes,
        "sync_rate": sync_successes / len(sync_records) if sync_records else math.nan,
        "timing_abs_mean": mean(timing_abs),
        "timing_abs_p95": percentile(timing_abs, 95),
        "acquisition_score_mean": mean(acq_scores),
        "runtime_ms_mean": mean(runtimes),
        "runtime_ms_p50": percentile(runtimes, 50),
        "runtime_ms_p95": percentile(runtimes, 95),
    }


def group_records(records: list[dict[str, Any]], keys: tuple[str, ...]) -> dict[tuple[Any, ...], list[dict[str, Any]]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for record in records:
        values = []
        for key in keys:
            if key == "profile":
                values.append(record["experiment"]["profile"])
            elif key == "payload_bits":
                values.append(record["signal"]["payload_bits"])
            elif key == "snr_db":
                values.append(record["channel"]["snr_db"])
            else:
                raise KeyError(key)
        grouped.setdefault(tuple(values), []).append(record)
    return grouped


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not math.isnan(float(value))


def fmt(value: Any, nd: int = 4) -> str:
    if value is None:
        return "nan"
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    if isinstance(value, (int, float)):
        return f"{float(value):.{nd}f}"
    return str(value)


def snr_threshold(points: list[tuple[float, float]], target: float) -> float:
    pts = sorted(points, key=lambda item: item[0])
    if not pts:
        return math.nan
    for idx in range(1, len(pts)):
        s0, r0 = pts[idx - 1]
        s1, r1 = pts[idx]
        if (r0 <= target <= r1) or (r1 <= target <= r0):
            if r1 == r0:
                return s1
            return s0 + (target - r0) * (s1 - s0) / (r1 - r0)
    if all(rate >= target for _, rate in pts):
        return pts[0][0]
    return math.nan


def write_markdown(path: Path, manifest: dict[str, Any], records: list[dict[str, Any]], summary_rows: list[dict[str, Any]]) -> None:
    global_summary = summarize(records)
    profile_rows = []
    for (profile,), recs in sorted(group_records(records, ("profile",)).items()):
        row = {"profile": profile}
        row.update(summarize(recs))
        profile_rows.append(row)

    threshold_rows = []
    grouped = group_records(records, ("profile", "payload_bits", "snr_db"))
    by_profile_payload: dict[tuple[str, int], list[tuple[float, float]]] = {}
    for (profile, payload_bits, snr_db), recs in grouped.items():
        by_profile_payload.setdefault((profile, payload_bits), []).append((float(snr_db), summarize(recs)["crc_rate"]))
    for (profile, payload_bits), points in sorted(by_profile_payload.items()):
        threshold_rows.append({
            "profile": profile,
            "payload_bits": payload_bits,
            "snr_at_50": snr_threshold(points, 0.50),
            "snr_at_90": snr_threshold(points, 0.90),
            "snr_at_95": snr_threshold(points, 0.95),
            "snr_at_99": snr_threshold(points, 0.99),
        })

    lines = [
        "# CDSS Simulate Monte Carlo Report",
        "",
        "## Campaign",
        "",
        f"- Created UTC: `{manifest['created_utc']}`",
        f"- Binary: `{manifest['binary']}`",
        f"- Binary SHA-256: `{manifest['binary_sha256']}`",
        f"- Total frames: **{global_summary['frames']}**",
        f"- Profiles: **{manifest['profiles']}**",
        f"- Payload bits: **{manifest['payloads']}**",
        f"- SNR points: **{manifest['snrs']}**",
        f"- Trials per point: **{manifest['trials_per_point']}**",
        "",
        "## Global Result",
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
    for row in profile_rows:
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
        "## Interpretation Boundaries",
        "",
        "- This campaign validates the public `cdss simulate` TX/channel/RX path.",
        "- Profiles using `--sync-search` also validate preamble acquisition and phase alignment.",
        "- Results are statistical evidence for the selected parameter grid, not a proof of universal robustness.",
        "- For publication, keep `manifest.json`, `results.jsonl`, CSV summaries, source revision, compiler, and binary hash together.",
        "",
        "## Output Files",
        "",
        "- `manifest.json`: campaign configuration and environment.",
        "- `results.jsonl`: one JSON object per simulated frame.",
        "- `summary_by_profile_payload_snr.csv`: main numeric performance table.",
        "- `summary_by_profile_payload.csv`: payload/profile aggregate table.",
        "- `summary_by_profile.csv`: profile aggregate table.",
        "- `failures.csv`: nonzero command exits or JSON parse failures.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_jobs(payloads: list[int], snrs: list[float], profiles: list[str], trials_per_point: int, seed_base: int) -> list[Job]:
    jobs: list[Job] = []
    job_id = 0
    for profile in profiles:
        if profile not in PROFILE_OPTIONS:
            raise ValueError(f"Unknown profile: {profile}; known profiles: {sorted(PROFILE_OPTIONS)}")
        for payload_bits in payloads:
            for snr_db in snrs:
                for trial in range(trials_per_point):
                    seed = seed_base + job_id
                    jobs.append(Job(job_id, profile, payload_bits, snr_db, trial, seed))
                    job_id += 1
    return jobs


def main() -> int:
    parser = argparse.ArgumentParser(description="CDSS Monte Carlo campaign using `cdss simulate`.")
    parser.add_argument("--preset", choices=sorted(PRESETS), default="smoke")
    parser.add_argument("--binary", help="Path or command name for the cdss executable")
    parser.add_argument("--out-dir", help="Output directory; default is timestamped under research/mc_results")
    parser.add_argument("--payloads", help="Comma-separated payload sizes in bits")
    parser.add_argument("--snrs", help="Comma-separated SNR points in dB")
    parser.add_argument("--profiles", help=f"Comma-separated profiles: {','.join(PROFILE_OPTIONS)}")
    parser.add_argument("--trials-per-point", type=int, help="Trials per profile/payload/SNR point")
    parser.add_argument("--seed-base", type=int, default=20260512)
    parser.add_argument("--workers", type=int, default=max(1, min(4, (os.cpu_count() or 2) // 2)))
    parser.add_argument("--omp-threads", type=int, default=1, help="OpenMP threads per cdss simulate process")
    parser.add_argument("--timeout-s", type=float, default=240.0)
    parser.add_argument("--fail-fast", action="store_true")
    args = parser.parse_args()

    preset = PRESETS[args.preset]
    payloads = parse_int_list(args.payloads) if args.payloads else list(preset["payloads"])
    snrs = parse_float_list(args.snrs) if args.snrs else list(preset["snrs"])
    profiles = parse_str_list(args.profiles) if args.profiles else list(preset["profiles"])
    trials_per_point = args.trials_per_point if args.trials_per_point is not None else int(preset["trials_per_point"])
    if trials_per_point <= 0:
        raise ValueError("--trials-per-point must be positive")

    binary = resolve_binary(args.binary)
    timestamp = datetime.now().strftime("simulate_mc_%Y%m%d_%H%M%S")
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (DEFAULT_OUT_ROOT / timestamp)
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(payloads, snrs, profiles, trials_per_point, args.seed_base)
    manifest = {
        "schema_version": "cdss-simulate-monte-carlo-manifest-v1",
        "created_utc": utc_now(),
        "preset": args.preset,
        "binary": str(binary),
        "binary_sha256": file_sha256(binary),
        "python": sys.version,
        "platform": platform.platform(),
        "repo_root": str(REPO_ROOT),
        "payloads": payloads,
        "snrs": snrs,
        "profiles": profiles,
        "profile_options": {name: PROFILE_OPTIONS[name] for name in profiles},
        "trials_per_point": trials_per_point,
        "seed_base": args.seed_base,
        "workers": args.workers,
        "omp_threads": args.omp_threads,
        "timeout_s": args.timeout_s,
        "total_jobs": len(jobs),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    results_path = out_dir / "results.jsonl"
    print(f"CDSS simulate Monte Carlo: {len(jobs)} frames")
    print(f"  binary: {binary}")
    print(f"  output: {out_dir}")

    records: list[dict[str, Any]] = []
    completed = 0
    failures = 0
    with results_path.open("w", encoding="utf-8") as f:
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
            future_map = {
                executor.submit(run_trial, binary, job, args.timeout_s, args.omp_threads): job
                for job in jobs
            }
            for future in as_completed(future_map):
                job = future_map[future]
                try:
                    record = future.result()
                except Exception as exc:  # noqa: BLE001 - stored as failed trial.
                    record = {
                        "schema_version": "cdss-simulate-monte-carlo-v1",
                        "experiment": {
                            "job_id": job.job_id,
                            "profile": job.profile,
                            "trial": job.trial,
                            "seed": job.seed,
                            "timestamp_utc": utc_now(),
                            "command": [],
                            "exit_code": -999,
                            "runtime_ms": math.nan,
                        },
                        "signal": {"payload_bits": job.payload_bits},
                        "channel": {"snr_db": job.snr_db},
                        "acquisition": {"sync_search": "--sync-search" in PROFILE_OPTIONS[job.profile], "sync_ok": False},
                        "decoder": {"decode_attempted": False, "crc_ok": False, "refinement_iters_used": 0, "payload_bit_count": 0},
                        "payload": {"bit_errors": job.payload_bits, "total_bits": job.payload_bits, "ber": 1.0},
                        "process": {"stdout": "", "stderr": str(exc), "parse_error": str(exc)},
                    }
                records.append(record)
                f.write(json.dumps(record, sort_keys=True) + "\n")
                f.flush()
                completed += 1
                if record["experiment"]["exit_code"] != 0 or record["process"].get("parse_error"):
                    failures += 1
                    if args.fail_fast:
                        raise RuntimeError(f"Trial failed: job_id={job.job_id}")
                if completed == len(jobs) or completed % max(1, min(25, len(jobs))) == 0:
                    crc_count = sum(1 for r in records if r["decoder"]["crc_ok"])
                    print(f"  completed {completed}/{len(jobs)}; crc_ok={crc_count}; failures={failures}", flush=True)

    summary_rows: list[dict[str, Any]] = []
    for key, recs in sorted(group_records(records, ("profile", "payload_bits", "snr_db")).items()):
        profile, payload_bits, snr_db = key
        row = {"profile": profile, "payload_bits": payload_bits, "snr_db": snr_db}
        row.update(summarize(recs))
        summary_rows.append(row)
    write_csv(out_dir / "summary_by_profile_payload_snr.csv", summary_rows)

    payload_rows: list[dict[str, Any]] = []
    for key, recs in sorted(group_records(records, ("profile", "payload_bits")).items()):
        profile, payload_bits = key
        row = {"profile": profile, "payload_bits": payload_bits}
        row.update(summarize(recs))
        payload_rows.append(row)
    write_csv(out_dir / "summary_by_profile_payload.csv", payload_rows)

    profile_rows: list[dict[str, Any]] = []
    for key, recs in sorted(group_records(records, ("profile",)).items()):
        (profile,) = key
        row = {"profile": profile}
        row.update(summarize(recs))
        profile_rows.append(row)
    write_csv(out_dir / "summary_by_profile.csv", profile_rows)

    failure_rows = []
    for record in records:
        if record["experiment"]["exit_code"] != 0 or record["process"].get("parse_error"):
            failure_rows.append({
                "job_id": record["experiment"]["job_id"],
                "profile": record["experiment"]["profile"],
                "seed": record["experiment"]["seed"],
                "payload_bits": record["signal"]["payload_bits"],
                "snr_db": record["channel"].get("snr_db"),
                "exit_code": record["experiment"]["exit_code"],
                "parse_error": record["process"].get("parse_error", ""),
                "stderr": record["process"].get("stderr", ""),
            })
    write_csv(out_dir / "failures.csv", failure_rows or [{"job_id": "", "profile": "", "seed": "", "payload_bits": "", "snr_db": "", "exit_code": "", "parse_error": "", "stderr": ""}])

    write_markdown(out_dir / "paper_summary.md", manifest, records, summary_rows)
    global_summary = summarize(records)
    print("Summary:")
    print(f"  CRC success: {global_summary['crc_successes']}/{global_summary['frames']} ({global_summary['crc_rate'] * 100:.2f}%)")
    print(f"  Aggregate BER: {global_summary['ber_aggregate']:.8f}")
    print(f"  Report: {out_dir / 'paper_summary.md'}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
