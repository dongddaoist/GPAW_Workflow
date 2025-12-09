#!/bin/bash
set -e

# BASE_DIR=$(pwd)
CSV="/fs/scratch/PYS1078/dd_scratch/materials-project/poscar-input-list.csv"
NUM=${1:-20}

# Select next NUM jobs and mark them as running
IDS=$(python3 select_next_jobs.py --base-dir "$BASE_DIR" --csv "$CSV" --num "$NUM")

if [ -z "$IDS" ]; then
  echo "No waiting jobs left."
  exit 0
fi

echo "Submitting the following material_ids:"
echo "$IDS"

while read -r MID; do
  [ -z "$MID" ] && continue
  echo "sbatch for $MID"
  sbatch --export=ALL,MID="$MID" gpaw_job.slurm
done <<< "$IDS"

