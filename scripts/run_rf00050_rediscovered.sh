#!/usr/bin/env bash
#SBATCH --partition=aoraki
#SBATCH --job-name=cie-crispr-candidate-multi
#SBATCH --output=cie_crispr_candidate_multi_%j.log
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate cie

set -euo pipefail

REPO_DIR="/projects/health_sciences/bms/biochemistry/gardner_group/daniel/cm-iterative-expansion"
CM="${REPO_DIR}/data/crispr-pos-control/candidates-and-rfam.cm"
FASTA="${REPO_DIR}/data/crispr-pos-control/CRISPR_Cas2_domain_flanking_downstream_filtered.fasta"
BINARY="${REPO_DIR}/bin/cmiterativeexpansion"
OUTPUT_DIR="${REPO_DIR}/results/drh_crispr_candidate_cm_expansion_${SLURM_JOB_ID}"

MAX_ROUNDS=10
MIN_NEW_HITS=20
EVALUE_THRESHOLD=1e-5
SEED=42

cd "$REPO_DIR"

for required_file in "$CM" "$FASTA"; do
    if [ ! -f "$required_file" ]; then
        echo "[slurm] ERROR: file not found: $required_file" >&2
        exit 1
    fi
done

if [ ! -f "$REPO_DIR/go.mod" ]; then
    echo "[slurm] ERROR: no go.mod found at $REPO_DIR" >&2
    exit 1
fi

echo "[slurm] job $SLURM_JOB_ID starting on $(hostname), $SLURM_CPUS_PER_TASK cpus, $(date)"
echo "[slurm] multi-CM input: $CM"
echo "[slurm] sequence pool:  $FASTA"
echo "[slurm] output:          $OUTPUT_DIR"

if [ ! -x "$BINARY" ]; then
    echo "[slurm] building $BINARY ..."
    mkdir -p "$(dirname "$BINARY")"
    go build -o "$BINARY" ./src/
fi

CMI_BINARY="$BINARY" scripts/run_multi_cm.sh \
    "$CM" \
    "$FASTA" \
    "$OUTPUT_DIR" \
    -max-rounds "$MAX_ROUNDS" \
    -min-new-hits "$MIN_NEW_HITS" \
    -evalue "$EVALUE_THRESHOLD" \
    -seed "$SEED" \
    -cpu "$SLURM_CPUS_PER_TASK"

echo "[slurm] done, $(date). Results in $OUTPUT_DIR"
