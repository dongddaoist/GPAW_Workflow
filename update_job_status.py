#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path
import os
import tempfile


def update_job_status(csv_path: Path, mid: str | None, status: str | None):
    """
    Ensure the CSV has a 'job_status' column.
    If mid and status are given, update that row's job_status.
    Writes back atomically via a temporary file + rename.
    """
    if not csv_path.exists():
        print(f"[ERR] CSV file not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    # Read entire CSV
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    # Ensure job_status column exists and is initialized
    if "job_status" not in fieldnames:
        fieldnames.append("job_status")
        for r in rows:
            # initialize only if missing
            if "job_status" not in r or r["job_status"] == "":
                r["job_status"] = "none"

    # Optionally update one row
    if mid is not None and status is not None:
        found = False
        for r in rows:
            if r.get("material_id", "").strip() == mid:
                r["job_status"] = status
                found = True
                break
        if not found:
            print(f"[WARN] material_id {mid} not found in {csv_path}", file=sys.stderr)

    # Write atomically
    tmp_path = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    os.replace(tmp_path, csv_path)


def main():
    p = argparse.ArgumentParser(
        description="Initialize or update job_status column in poscar-input-list.csv"
    )
    p.add_argument("--csv", required=True, help="Path to poscar-input-list.csv")
    p.add_argument("--mid", help="material_id whose status to update")
    p.add_argument("--status", help="New job_status value (e.g., none/waiting/running/completed)")

    args = p.parse_args()
    csv_path = Path(args.csv).resolve()

    # If mid or status missing, we just ensure the column exists and exit.
    mid = args.mid
    status = args.status
    if (mid is None) != (status is None):
        print("[ERR] --mid and --status must be given together, or both omitted.", file=sys.stderr)
        sys.exit(1)

    update_job_status(csv_path, mid, status)


if __name__ == "__main__":
    main()

