#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import sys
import os

def load_table(csv_path: Path):
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return rows, fieldnames


def save_table_atomic(csv_path: Path, rows, fieldnames):
    tmp_path = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})
    os.replace(tmp_path, csv_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", default=".", help="project root (contains poscar-input/)")
    ap.add_argument("--csv", default="poscar-input/poscar-input-list.csv")
    ap.add_argument("--num", type=int, default=20,
                    help="number of jobs to select")
    args = ap.parse_args()

    base = Path(args.base_dir).resolve()
    csv_path = base / args.csv
    poscar_input_dir = base / "poscar-input"

    if not csv_path.exists():
        print(f"[ERR] CSV {csv_path} not found", file=sys.stderr)
        sys.exit(1)

    rows, fieldnames = load_table(csv_path)
    if "job_status" not in fieldnames:
        fieldnames.append("job_status")

    # Initialize missing job_status entries to 'none'
    for r in rows:
        if r.get("job_status", "").strip() == "":
            r["job_status"] = "none"

    # 1) refresh job_status from DONE markers
    for r in rows:
        mid = r.get("material_id", "").strip()
        if not mid:
            continue
        mp_dir = poscar_input_dir / mid
        done_file = mp_dir / "DONE"
        if done_file.exists():
            r["job_status"] = "completed"

    # 2) select first N that are not completed or running
    selected = []
    count = 0
    for r in rows:
        if count >= args.num:
            break
        status = r.get("job_status", "").strip().lower()
        if status in ("completed", "running"):
            continue
        mid = r.get("material_id", "").strip()
        if not mid:
            continue
        # mark as running
        r["job_status"] = "running"
        selected.append(mid)
        count += 1

    # 3) write back updated CSV atomically
    save_table_atomic(csv_path, rows, fieldnames)

    # 4) print selected IDs for bash to consume (one per line)
    for mid in selected:
        print(mid)


if __name__ == "__main__":
    main()

