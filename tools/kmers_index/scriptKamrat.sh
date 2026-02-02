#!/usr/bin/env bash
set -euo pipefail

TIME="${TIME:-}"

KAMRAT_IMG="$1"
IN_MAT="$2"
OUT_DIR="$3"
kmerSize="$4"

KAMRAT_NORMALIZE="${5:-false}"   # true/false
NFBASE="${6:-1000000}"

IDX_DIR="$OUT_DIR/index"
mkdir -p "$IDX_DIR" "$OUT_DIR"

kmerMaxOverlap=$((kmerSize - 1))
kmerMinOverlap=$((kmerSize / 2))

# Optional args for kamrat index
INDEX_ARGS=()

# normalization
if [[ "${KAMRAT_NORMALIZE,,}" == "true" ]]; then
  INDEX_ARGS+=(-nfbase "$NFBASE")
fi

$TIME apptainer exec -B '/store:/store' -B '/data:/data' "$KAMRAT_IMG" \
  kamrat index \
    -intab "$IN_MAT" \
    -outdir "$IDX_DIR" \
    -klen "$kmerSize" \
    "${INDEX_ARGS[@]}"

$TIME apptainer exec -B '/store:/store' -B '/data:/data' "$KAMRAT_IMG" \
  kamrat merge \
    -idxdir "$IDX_DIR" \
    -overlap "${kmerMaxOverlap}-${kmerMinOverlap}" \
    -interv pearson:0.20 \
    -counts mean:float \
    -outfmt tab \
    -outpath "$OUT_DIR/merged-res.tsv"
