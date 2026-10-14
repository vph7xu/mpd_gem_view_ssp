#!/usr/bin/env python3
"""Batch CERN GEM viewer replay over both UVA_DATA trees.

Usage: python3 analyze_all_runs.py --write-calibrations run_calibrations.csv
       python3 analyze_all_runs.py --calibrations run_calibrations.csv --slurm --dry-run
       python3 analyze_all_runs.py --calibrations run_calibrations.csv --slurm --account YOUR_ALLOCATION
Run from the mpd_gem_view_ssp repository root. Paths in the CSV are relative
to that root; EVIO paths are discovered recursively below --data-root.
For Rivanna, pass --setup setup_env.sh if your replay requires environment setup.
Each job runs all files of one run in sequence; --plot adds one ROOT plot at the end.
Outputs are collected in Rootfiles/<calibrations CSV filename without extension>/.
For example, prefix.csv gives Rootfiles/prefix/ for every selected run.
"""
import argparse
import csv
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def run_number(path):
    # Accept run472, run_472, and exp1_472; avoid matching dates, HV and segments.
    name = path.name
    matches = re.findall(r"(?i)(?:^|[^a-z0-9])run[_-]?(\d{3,5})(?!\d)", name)
    if not matches:
        matches = re.findall(r"(?i)(?:^|[^a-z0-9])exp\d+[_-](\d{3,5})(?!\d)", name)
    if len(set(matches)) != 1:
        return None
    return int(matches[0])


def segment_number(path):
    match = re.search(r"\.evio\.(\d+)$", path.name, re.I)
    return int(match.group(1)) if match else -1


def read_calibrations(path, repo):
    result = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(row for row in stream if row.strip() and not row.lstrip().startswith("#"))
        expected = {"run", "pedestal", "common_mode"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError("CSV needs headers: run,pedestal,common_mode")
        for line, row in enumerate(reader, 2):
            spec = (row["run"] or "").strip()
            match = re.fullmatch(r"(\d+)(?:\s*-\s*(\d+))?", spec)
            if not match:
                raise ValueError(f"CSV row {line}: invalid run/range {spec!r}")
            first = int(match.group(1))
            last = int(match.group(2) or first)
            if last < first or last - first > 10000:
                raise ValueError(f"CSV row {line}: invalid range {spec!r}")
            ped = (row["pedestal"] or "").strip()
            cm = (row["common_mode"] or "").strip()
            if not ped or not cm:
                raise ValueError(f"CSV row {line}: missing calibration path")
            paths = tuple((repo / p).resolve() if not Path(p).is_absolute() else Path(p) for p in (ped, cm))
            for number in range(first, last + 1):
                if number in result:
                    raise ValueError(f"CSV row {line}: duplicate calibration for run {number}")
                result[number] = paths
    return result


def replay_command(replay, ped, cm, path):
    return [str(replay), "-c", "0", "-t", "0", "-z", "1", "-n", "-1", "--tracking", "off",
            "--pedestal", str(ped), "--common_mode", str(cm), str(path)]


def calibration_output_dir(repo, calibration_csv):
    prefix = Path(calibration_csv).stem
    if prefix in ("", ".", ".."):
        raise ValueError("invalid calibrations CSV filename for output directory")
    return repo / "Rootfiles" / prefix


def collect_run_outputs(repo, number, calibration_csv):
    """Move this run's files only; other concurrent run outputs remain untouched."""
    source = repo / "Rootfiles"
    destination = calibration_output_dir(repo, calibration_csv)
    destination.mkdir(parents=True, exist_ok=True)
    count = 0
    for output in source.iterdir():
        if output.is_file() and run_number(output) == number:
            shutil.move(str(output), str(destination / output.name))
            count += 1
    print(f"Run {number}: moved {count} output file(s) to {destination}", flush=True)
    if count == 0:
        print(f"Warning: no output files with run {number} in their filename were found in {source}.", file=sys.stderr)
    return destination


def plot_command(repo, number, calibration_csv):
    output = calibration_output_dir(repo, calibration_csv) / f"cluster_0_cern_beamtest_exp1_{number}.root_data_quality_check.root"
    macro = f'scripts/plot_quality_results_2.cpp("{output}")'
    return output, ["root", "-b", "-q", macro]


def write_slurm_job(path, repo, replay, number, files, ped, cm, args):
    lines = ["#!/bin/bash", f"#SBATCH --job-name=gem_{number}",
             f"#SBATCH --partition={args.partition}",
             f"#SBATCH --time={args.time}",
             f"#SBATCH --cpus-per-task={args.cpus}",
             f"#SBATCH --mem={args.mem}",
             f"#SBATCH --output={args.log_dir.resolve()}/slurm_run{number}_%j.out"]
    if args.account:
        lines.append(f"#SBATCH --account={args.account}")
    lines += ["", "set -euo pipefail", f"cd {shlex.quote(str(repo))}"]
    if args.setup:
        lines.append(f"source {shlex.quote(str(args.setup.resolve()))}")
    lines.append(f"export OMP_NUM_THREADS={args.cpus}")
    for raw in files:
        lines.append(shlex.join(replay_command(replay, ped, cm, raw)))
    lines.append(shlex.join(["python3", str(Path(__file__).resolve()), "--repo", str(repo),
                            "--collect-output", str(number), str(args.calibrations.resolve())]))
    if args.plot:
        output, cmd = plot_command(repo, number, args.calibrations)
        lines += [f"test -f {shlex.quote(str(output))}", shlex.join(cmd)]
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, default=Path("/standard/sbs/CERN_BEAM_TEST"))
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--collect-output", nargs=2, metavar=("RUN", "CALIBRATION_CSV"), help=argparse.SUPPRESS)
    parser.add_argument("--replay", type=Path, help="replay executable; defaults to repo/bin/replay or repo/replay")
    parser.add_argument("--calibrations", type=Path, help="run to calibration CSV")
    parser.add_argument("--write-calibrations", type=Path, metavar="CSV",
                        help="create a blank row for each discovered run, then exit; never overwrite")
    parser.add_argument("--runs", help="optional runs, e.g. 293,294,472-481")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--plot", action="store_true", help="run ROOT quality macro after each successful replay")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--log-dir", type=Path, default=Path("batch_logs"))
    parser.add_argument("--slurm", action="store_true", help="create and submit one sbatch job per run")
    parser.add_argument("--account", help="Slurm allocation/account, if required")
    parser.add_argument("--partition", default="standard")
    parser.add_argument("--time", default="12:00:00", help="wall time per run")
    parser.add_argument("--mem", default="16G", help="memory per run")
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--setup", type=Path, help="optional environment setup script sourced in each job")
    args = parser.parse_args()
    if args.collect_output:
        number, calibration_csv = args.collect_output
        collect_run_outputs(args.repo.resolve(), int(number), calibration_csv)
        return 0
    if not args.calibrations and not args.write_calibrations:
        parser.error("specify --calibrations or --write-calibrations")
    if args.cpus < 1 or not re.fullmatch(r"[A-Za-z0-9_-]+", args.partition):
        parser.error("invalid --cpus or --partition")
    if args.account and not re.fullmatch(r"[A-Za-z0-9_-]+", args.account):
        parser.error("invalid --account")
    if not re.fullmatch(r"(?:\d+-)?\d{1,2}:\d{2}:\d{2}", args.time):
        parser.error("invalid --time; use HH:MM:SS or D-HH:MM:SS")
    if not re.fullmatch(r"\d+[KMGTP]?", args.mem, re.I):
        parser.error("invalid --mem; example: 8G")
    if args.setup and not args.setup.resolve().is_file():
        parser.error(f"setup script missing: {args.setup}")
    repo = args.repo.resolve()
    data_root = args.data_root.resolve()
    selected = None
    if args.runs:
        selected = set()
        for token in args.runs.split(","):
            match = re.fullmatch(r"\s*(\d+)(?:\s*-\s*(\d+))?\s*", token)
            if not match:
                parser.error(f"invalid --runs token: {token!r}")
            selected.update(range(int(match.group(1)), int(match.group(2) or match.group(1)) + 1))

    groups = defaultdict(list)
    unknown = []
    for directory in (data_root / "UVA_DATA_1", data_root / "UVA_DATA_2"):
        if not directory.is_dir():
            parser.error(f"missing data directory: {directory}")
        for path in directory.rglob("*"):
            if not path.is_file() or not re.search(r"\.evio\.\d+$", path.name, re.I):
                continue
            number = run_number(path)
            if number is None:
                unknown.append(path)
            elif selected is None or number in selected:
                groups[number].append(path)
    if unknown:
        print(f"Skipped {len(unknown)} EVIO files without an unambiguous run number; first examples:", file=sys.stderr)
        for path in unknown[:10]:
            print(f"  {path}", file=sys.stderr)
    if not groups:
        parser.error("no matching run EVIO files found")
    if args.write_calibrations:
        destination = args.write_calibrations
        try:
            with destination.open("x", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("run", "pedestal", "common_mode"))
                for number in sorted(groups):
                    writer.writerow((number, "", ""))
        except FileExistsError:
            parser.error(f"CSV already exists, refusing to overwrite: {destination}")
        print(f"Wrote {len(groups)} run rows to {destination}. Fill pedestal and common_mode before replay.")
        return 0
    calibration = read_calibrations(args.calibrations, repo)
    # The CSV is the allow-list. Unlisted EVIO runs are never replayed or submitted.
    unlisted = set(groups) - set(calibration)
    absent = set(calibration) - set(groups)
    if selected is not None:
        absent &= selected
    groups = {number: files for number, files in groups.items() if number in calibration}
    print(f"Selected {len(groups)} run(s) present in CSV and data directories; "
          f"ignored {len(unlisted)} unlisted data run(s).")
    if absent:
        print(f"CSV run(s) without EVIO files ({len(absent)}): " + ", ".join(map(str, sorted(absent))),
              file=sys.stderr)
    if not groups:
        print("No CSV-listed runs have EVIO files. Nothing to submit.", file=sys.stderr)
        return 0
    for number in groups:
        for path in calibration[number]:
            if not path.is_file():
                parser.error(f"run {number}: calibration file missing: {path}")
    if args.replay:
        replay = (repo / args.replay).resolve()
    else:
        replay = next((p for p in (repo / "bin/replay", repo / "replay") if p.is_file()), repo / "bin/replay")
    if not replay.is_file():
        parser.error(f"replay binary missing: {replay}; specify its location with --replay")
    if args.plot and not (repo / "scripts/plot_quality_results_2.cpp").is_file():
        parser.error("scripts/plot_quality_results_2.cpp missing")

    # Replay writes to Rootfiles by input basename. Reject basename collisions
    # before processing so one raw file cannot overwrite another's output.
    names = defaultdict(list)
    for files in groups.values():
        for path in files:
            names[path.name].append(path)
    collisions = {name: paths for name, paths in names.items() if len(paths) > 1}
    if collisions:
        parser.error("duplicate EVIO basenames across data trees: " + ", ".join(sorted(collisions)[:10]))

    if not args.dry_run or args.slurm:
        args.log_dir.mkdir(parents=True, exist_ok=True)
    failures = []
    for number in sorted(groups):
        ped, cm = calibration[number]
        files = sorted(groups[number], key=lambda p: (p.name.rsplit(".evio.", 1)[0], segment_number(p)))
        print(f"\nRUN {number}: {len(files)} file(s), pedestal={ped}, common_mode={cm}", flush=True)
        print(f"Output directory: {calibration_output_dir(repo, args.calibrations)}", flush=True)
        if args.slurm:
            job = args.log_dir.resolve() / f"run{number}.sbatch"
            write_slurm_job(job, repo, replay, number, files, ped, cm, args)
            command = ["sbatch", "--parsable", str(job)]
            print(f"Job script: {job}\n{shlex.join(command)}", flush=True)
            if not args.dry_run:
                result = subprocess.run(command, cwd=repo, text=True, capture_output=True)
                if result.returncode:
                    failures.append((str(job), result.returncode, result.stderr.strip()))
                    print(f"SUBMISSION FAILED: {result.stderr.strip()}", file=sys.stderr)
                    if not args.continue_on_error:
                        return 1
                else:
                    print(f"Submitted run {number}: job {result.stdout.strip()}", flush=True)
            continue
        for path in files:
            cmd = replay_command(replay, ped, cm, path)
            print(shlex.join(cmd), flush=True)
            if args.dry_run:
                continue
            log = args.log_dir / f"run{number}_{path.name}.log"
            with log.open("w") as stream:
                status = subprocess.run(cmd, cwd=repo, stdout=stream, stderr=subprocess.STDOUT).returncode
            if status:
                failures.append((str(path), status, str(log)))
                print(f"FAILED (exit {status}): {log}", file=sys.stderr)
                if not args.continue_on_error:
                    return 1
                continue
        if not args.dry_run:
            collect_run_outputs(repo, number, args.calibrations)
        if args.plot and not args.dry_run:
            output, plot_cmd = plot_command(repo, number, args.calibrations)
            if not output.is_file():
                failures.append((str(output), 1, "ROOT output missing"))
                print(f"Plot failed: ROOT output missing: {output}", file=sys.stderr)
                if not args.continue_on_error:
                    return 1
            else:
                plot_log = args.log_dir / f"run{number}.plot.log"
                with plot_log.open("w") as stream:
                    status = subprocess.run(plot_cmd, cwd=repo,
                                            stdout=stream, stderr=subprocess.STDOUT).returncode
                if status:
                    failures.append((str(output), status, str(plot_log)))
                    if not args.continue_on_error:
                        return 1
    print(f"\nCompleted {len(groups)} run(s); {len(failures)} failure(s).")
    for name, status, log in failures:
        print(f"  exit {status}: {name}  log: {log}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
