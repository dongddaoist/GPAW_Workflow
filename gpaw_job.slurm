#!/bin/bash
#SBATCH --time=96:00:00 # walltime, abbreviated by -t
#SBATCH --ntasks=40  # number of MPI tasks, abbreviated by -n
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --account=PYS1137 # account - abbreviated by -A
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=cherrydons@gmail.com
#SBATCH -o mdex.out
#SBATCH -e mdex.err

JOB_START_TIME=$(date +%s)


# Make sure logs/ exists
mkdir -p logs

if [ -z "$MID" ]; then
  echo "Error: MID not set (material_id)" >&2
  exit 1
fi

BASE_DIR=$(pwd)

module load miniconda3/24.1.2-py310
eval "$(conda shell.bash hook)"
conda activate gpaw

# Use MPI parallelization
srun -n 40 python3 gpaw_workflow_single.py \
    --mid "$MID" \
    --base-dir "$BASE_DIR" \
    --mode fd \
    --h 0.18 \
    --gridref 4 \
#    --compute_elastic

JOB_END_TIME=$(date +%s)
TOTAL_SECONDS=$((JOB_END_TIME - JOB_START_TIME))

# -------------------------------
# RECORD TIMING (SLURM-LEVEL)
# -------------------------------
TRACKER_FILE="${BASE_DIR}/time_tracker.txt"

{
  echo "${MID}, ${TOTAL_SECONDS} seconds, jobid=${SLURM_JOB_ID}, ntasks=${SLURM_NTASKS}"
} >> "$TRACKER_FILE"

echo "[SLURM TIME] ${MID} completed in ${TOTAL_SECONDS} seconds"

