#!/usr/bin/env bash
# Run cmiterativeexpansion independently for every model in a multi-CM file.
# Usage: run_multi_cm.sh models.cm sequences.fasta output_dir [cmiterativeexpansion flags...]

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: run_multi_cm.sh MULTI_CM FASTA OUTPUT_DIR [cmiterativeexpansion flags...]

The first three arguments are consumed by this wrapper. Remaining arguments are
passed to cmiterativeexpansion, so do not pass -cm, -fasta, or -out there.

Environment:
  CMI_BINARY  cmiterativeexpansion executable (default: cmiterativeexpansion)
EOF
    exit 2
}

if [ "$#" -lt 3 ]; then
    usage
fi

CM_FILE=$1
FASTA=$2
OUTPUT_DIR=$3
shift 3

CMI_BINARY=${CMI_BINARY:-cmiterativeexpansion}

for required_file in "$CM_FILE" "$FASTA"; do
    if [ ! -f "$required_file" ]; then
        echo "[multi-cm] ERROR: file not found: $required_file" >&2
        exit 1
    fi
done

for required_command in cmstat cmfetch "$CMI_BINARY"; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "[multi-cm] ERROR: command not found: $required_command" >&2
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"
MANIFEST="$OUTPUT_DIR/models.tsv"
WORK_DIR="$OUTPUT_DIR/.extracted_models"
mkdir -p "$WORK_DIR"

printf 'index\tmodel_name\tmodel_file\toutput_dir\n' > "$MANIFEST"

model_count=0
while IFS= read -r model_name; do
    [ -n "$model_name" ] || continue
    model_count=$((model_count + 1))

    # Prefixing with the model index prevents collisions after sanitization.
    safe_name=$(printf '%s' "$model_name" | tr -c '[:alnum:]_.-' '_')
    model_id=$(printf 'model_%03d_%s' "$model_count" "$safe_name")
    model_file="$WORK_DIR/$model_id.cm"
    model_output="$OUTPUT_DIR/$model_id"

    echo "[multi-cm] [$model_count] extracting $model_name"
    if ! cmfetch "$CM_FILE" "$model_name" > "$model_file"; then
        echo "[multi-cm] ERROR: cmfetch failed for model $model_name" >&2
        exit 1
    fi
    if [ ! -s "$model_file" ]; then
        echo "[multi-cm] ERROR: cmfetch returned an empty model for $model_name" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\t%s\n' \
        "$model_count" "$model_name" "$model_file" "$model_output" >> "$MANIFEST"

    echo "[multi-cm] [$model_count] expanding $model_name"
    "$CMI_BINARY" -cm "$model_file" -fasta "$FASTA" -out "$model_output" "$@"
done < <(cmstat "$CM_FILE" | awk '$1 ~ /^[0-9]+$/ { print $2 }')

if [ "$model_count" -eq 0 ]; then
    echo "[multi-cm] ERROR: cmstat found no models in $CM_FILE" >&2
    exit 1
fi

echo "[multi-cm] completed $model_count models"
echo "[multi-cm] manifest: $MANIFEST"
