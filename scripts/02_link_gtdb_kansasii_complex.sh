#!/usr/bin/env bash

# srun --pty -n 1 -c 6 --time=01:00:00 --mem=16G bash -l
#
# Create species-grouped symlinks to the already-downloaded genomes of the
# M. kansasii species complex, based on the accessions list produced by
# scripts/01_download_filter_gtdb.sh.
#
# Usage: 02_link_gtdb_kansasii_complex.sh [--accessions FILE] [--ncbi-dataset DIR] [--link-root DIR]
#   defaults: GTDB_ACCESSIONS_DIR/kansasii-complex-mlsa-accessions.txt, NCBI_DATASET, GTDB_COMPLEX_GENOMES

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"

ACCESSIONS_OUT=$GTDB_ACCESSIONS_DIR/kansasii-complex-mlsa-accessions.txt
NCBI_DATA=$NCBI_DATASET
LINK_ROOT=$GTDB_COMPLEX_GENOMES
while [ $# -gt 0 ]; do
  case "$1" in
    --accessions) ACCESSIONS_OUT=$2; shift 2 ;;
    --ncbi-dataset) NCBI_DATA=$2; shift 2 ;;
    --link-root) LINK_ROOT=$2; shift 2 ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac
done

if [ ! -f "$ACCESSIONS_OUT" ]; then
  echo "ERROR: $ACCESSIONS_OUT not found; run scripts/01_download_filter_gtdb.sh first" >&2
  exit 1
fi

# --- create species-grouped symlinks to the already-downloaded genomes -------

mkdir -p "$LINK_ROOT"
for species in kansasii persicum pseudokansasii innocens attenuatum ostraviense gastri; do
  mkdir -p "$LINK_ROOT/$species"
done

n_linked=0
n_missing=0
while IFS=$'\t' read -r acc species; do
  gdir="$NCBI_DATA/$acc"
  if [ ! -d "$gdir" ]; then
    echo "WARNING: no local genome directory for $acc (species: $species), skipping" >&2
    n_missing=$((n_missing + 1))
    continue
  fi
  fna=$(find "$gdir" -maxdepth 1 -iname '*_genomic.fna' | head -1)
  gff=$(find "$gdir" -maxdepth 1 -iname 'genomic.gff' | head -1)
  if [ -z "$fna" ] || [ -z "$gff" ]; then
    echo "WARNING: missing .fna or genomic.gff for $acc (species: $species), skipping" >&2
    n_missing=$((n_missing + 1))
    continue
  fi
  ln -sfn "$fna" "$LINK_ROOT/$species/${acc}.fna"
  ln -sfn "$gff" "$LINK_ROOT/$species/${acc}.gff"
  n_linked=$((n_linked + 1))
done < "$ACCESSIONS_OUT"

echo "Linked $n_linked genomes ($n_missing missing/skipped) into $LINK_ROOT/<species>/"
echo "Per-species genome counts in $LINK_ROOT:"
for species in kansasii persicum pseudokansasii innocens attenuatum ostraviense gastri; do
  count=$(find "$LINK_ROOT/$species" -iname '*.fna' | wc -l)
  echo "  $species: $count"
done
