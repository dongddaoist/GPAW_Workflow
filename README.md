# GPAW_highthroughput
High-Throughput GPAW + Bader Charge Pipeline for Materials Project Structures
Overview

This repository provides a high-throughput, MPI-parallel workflow for:
Computing accurate all-electron Bader charges using GPAW
Handling both optimized and non-optimized structures
Optionally computing elastic tensors (6×6 Voigt matrix)
Running at scale on HPC systems with Slurm

Tracking job status, completion, and total Slurm wall time

The pipeline is designed for large materials datasets (O(10³–10⁴) structures), such as those retrieved from the Materials Project, and is suitable for force-field training, ML potentials, and charge models.

Key Features

✅ MPI-parallel GPAW (one structure per node) \
✅ Automatic handling of:

optimized geometries (single-point SCF)

non-optimized geometries (geometry relaxation + SCF)

✅ All-electron density (density_ae.cube) for reliable Bader analysis \
✅ Henkelman-group Bader charge partitioning \\ 
✅ Optional elastic tensor calculation \
✅ Robust job tracking (waiting → running → completed) \
✅ Automatic job submission and throttling \
✅ Total Slurm wall-time tracking per structure \

Directory Structure
project_root/
├── poscar-input/
│   ├── poscar-input-list.csv
│   └── mp-XXXXX/
│       └── POSCAR
├── retrieve-Li-O/           # optional (original MP downloads)
├── retrieve-Li-S/
│
├── gpaw_workflow_single.py  # main GPAW + Bader workflow
├── gpaw_job.slurm           # Slurm job script (1 structure / node)
├── submit_next_jobs.sh      # batch submission driver
├── select_next_jobs.py      # job selector / state updater
├── auto_submit.sh           # (optional) cron-based automation
│
├── logs/
├── time_tracker.txt         # appended Slurm wall time log
└── README.md

Input Data
1. POSCAR files

Each structure must exist as:
poscar-input/mp-XXXXXX/POSCAR

2. Master metadata file: poscar-input-list.csv

Required columns:

Column	Description
material_id	Materials Project ID (e.g. mp-1176651)
system	System label (e.g. Li-O, Li-S)
source_type	optimized or non_optimized
tensor	true / false (elastic tensor availability or request)

Optional / auto-added columns:

Column	Meaning
job_status	waiting, running, completed
Software Requirements
Core dependencies

Python ≥ 3.9

GPAW (MPI enabled)

ASE
mpi4py
NumPy
Slurm (HPC scheduler)
Bader charge analysis (REQUIRED)
The Henkelman-group bader executable must be available on compute nodes.

Supported options:
conda install -c conda-forge bader

OR manual compilation from:
https://theory.cm.utexas.edu/henkelman/code/bader/

Verify with:

bader --help
How the Workflow Works (Conceptual)

For each material (mp-XXXX):
Read POSCAR

If source_type == non_optimized:
Geometry relaxation (fixed cell)
Single-point SCF with GPAW
Compute all-electron charge density
Run Bader partitioning (rank-0 only)
Save per-atom Bader charges
(Optional) Compute elastic tensor
Mark job as completed
Log total Slurm wall time
Each material runs as an independent Slurm job using one node.
Running a Single Test Job (Strongly Recommended)
Before batch submission:
srun -n 1 python3 gpaw_workflow_single.py \
  --mid mp-1176651 \
  --base-dir /path/to/project_root \
  --mode fd \
  --h 0.18 \
  --gridref 4 \
  --compute_elastic


Check outputs under:

poscar-input/mp-1176651/gpaw/

Slurm Job Script

Each Slurm job:

Uses 1 node

Uses 40 MPI ranks

Sets OMP_NUM_THREADS=1

Runs one material ID via environment variable MID

Example submission:

MID=mp-1176651 sbatch gpaw_job.slurm

Batch Submission (High-Throughput Mode)

Submit the next batch (default: 20):

./submit_next_jobs.sh


Submit a smaller batch:

./submit_next_jobs.sh 5


The script will:

Detect completed jobs via DONE markers

Select first N waiting jobs

Mark them running

Submit one Slurm job per material

Automated Scheduling (Optional)

For continuous operation, use a cron job:

*/15 * * * * /bin/bash /path/to/project_root/auto_submit.sh


Includes a lock file to prevent concurrent submissions.

Outputs per Material

Located in:

poscar-input/mp-XXXX/gpaw/

File	Description
calc.gpw	GPAW restart file
POSCAR_final	Final geometry
density_ae.cube	All-electron density
density_valence.cube	Valence density (if available)
ACF.dat	Raw Bader output
bader_charges.csv	Per-atom charges
bader_summary.json	Full structured results
elastic_tensor_gpaw.txt	6×6 elastic tensor (optional)

Completion markers:

poscar-input/mp-XXXX/DONE

Time Tracking

Each job appends one line to:

time_tracker.txt


Format:

mp-1176651, 36821 seconds, jobid=42400404, ntasks=40


This reflects true Slurm wall time including:

environment loading

data I/O

GPAW computation

Bader analysis

MPI Safety Notes (Important)

GPAW is MPI-parallel → all ranks participate

External programs (bader, file writes):

Only rank 0 executes

Global synchronization uses:

world.barrier()


Failure to respect this will cause race conditions or MPI aborts.

Failure Recovery

If a job crashes:

Inspect Slurm logs in logs/

Inspect gpaw.log

Fix the issue

Reset status manually:

Remove RUNNING

Set job_status=waiting in CSV

Resubmit

Recommended Use Cases

Force-field training (Bader reference charges)

Charge-equilibration ML models

Benchmarking charge schemes

Elastic-property datasets

High-throughput materials screening

License

Specify as appropriate (e.g. MIT, BSD, internal use).

Acknowledgments

GPAW developers

ASE developers

Henkelman group (Bader analysis)

Materials Project
