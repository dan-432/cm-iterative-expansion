#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=cmiterativeexpansion
#SBATCH --output=cmsimile_%j.log
#SBATCH --time=35:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate cie

set -euo pipefail
# -------------------------------------------------------------------------
# Config -- edit these per run
# -------------------------------------------------------------------------

REPO_DIR="/projects/health_sciences/bms/biochemistry/gardner_group/daniel/cm-iterative-expansion"
CM="${REPO_DIR}/data/RF00050_1.cacofold.cm"
FASTA="${REPO_DIR}/data/discovery_pipeline/cluster_1341_domain_flanking_upstream_filtered.fasta"
BINARY="${REPO_DIR}/bin/cmiterativeexpansion"
BENCH_DIR="${REPO_DIR}/results/drh_fmntrue_cm_expansion_${SLURM_JOB_ID}"

MAX_ITERATIONS=10
MIN_NEW_HITS=20
EVALUE_THRESHOLD=1e-5
SEED=42
# -------------------------------------------------------------------------
# Environment
# -------------------------------------------------------------------------

cd "$REPO_DIR"

# Fail fast with a clear message rather than a cryptic Go toolchain error
# if REPO_DIR is ever wrong again (repo moved, script copied elsewhere,
# submitted from an unexpected location, etc.).
if [ ! -f "$REPO_DIR/go.mod" ]; then
    echo "[slurm] ERROR: no go.mod found at $REPO_DIR -- REPO_DIR is wrong. Check the hardcoded path at the top of this script." >&2
    exit 1
fi

echo "[slurm] job $SLURM_JOB_ID starting on $(hostname), $SLURM_CPUS_PER_TASK cpus, $(date)"
# -------------------------------------------------------------------------
# Build the binary if it isn't already there. Cheap and idempotent, so
# it's safe to leave in even once bin/ is committed/cached between runs.
# -------------------------------------------------------------------------
if [ ! -x "$BINARY" ]; then
    echo "[slurm] building $BINARY ..."
    mkdir -p "$(dirname "$BINARY")"
    go build -o "$BINARY" src/cmiterativeexpansion.go
fi

mkdir -p "$BENCH_DIR/combined_models"


# -------------------------------------------------------------------------
# run the cmiterativeexpansion with combination mode over that population.
# -------------------------------------------------------------------------

"$BINARY" -cm "$CM" \
    -fasta "$FASTA" -max-rounds "$MAX_ITERATIONS" \
    -min-new-hits "$MIN_NEW_HITS" -evalue "$EVALUE_THRESHOLD" \
    -seed "$SEED" -cpu "$SLURM_CPUS_PER_TASK" \
    -out "$BENCH_DIR"


echo "[slurm] done, $(date). Results in $BENCH_DIR"