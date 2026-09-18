#!/usr/bin/env bash

# srun --pty -n 1 -c 6 --time=01:00:00 --mem=16G bash -l

set -euo pipefail
mkdir -p /shares/sander.imm.uzh/MM/kansasii/lit/gtdb/gtdb232
cd /shares/sander.imm.uzh/MM/kansasii/lit/gtdb/gtdb232

if [ ! -f bac120_metadata_r232.tsv ]; then
  wget https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/bac120_metadata_r232.tsv.gz
  gunzip bac120_metadata_r232.tsv.gz
fi

# --- filter GTDB metadata down to the M. kansasii species complex ------------

SPECIES_PATTERN='s__Mycobacterium (kansasii|persicum|pseudokansasii|innocens|attenuatum|ostraviense|gastri)'
METADATA=bac120_metadata_r232.tsv
ACCESSIONS_OUT=kansasii-complex-mlsa-accessions.txt

# Column indices are looked up from the header at runtime so this stays
# robust to GTDB reordering columns between releases.
header=$(head -1 "$METADATA")
acc_col=$(echo "$header" | tr '\t' '\n' | grep -nx 'accession' | cut -d: -f1)
tax_col=$(echo "$header" | tr '\t' '\n' | grep -nx 'gtdb_taxonomy' | cut -d: -f1)

if [ -z "$acc_col" ] || [ -z "$tax_col" ]; then
  echo "ERROR: could not find 'accession' or 'gtdb_taxonomy' columns in $METADATA" >&2
  exit 1
fi

awk -F'\t' -v acc_col="$acc_col" -v tax_col="$tax_col" -v pattern="$SPECIES_PATTERN" '
  NR > 1 && $tax_col ~ pattern {
    acc = $acc_col
    sub(/^(RS_|GB_)/, "", acc)
    tax = $tax_col
    sub(/.*s__Mycobacterium /, "", tax)
    print acc "\t" tax
  }
' "$METADATA" | sort -u > "$ACCESSIONS_OUT"

n_acc=$(wc -l < "$ACCESSIONS_OUT")
echo "Wrote $n_acc kansasii-complex accessions to $PWD/$ACCESSIONS_OUT"
echo "Per-species counts:"
cut -f2 "$ACCESSIONS_OUT" | sort | uniq -c

# --- create species-grouped symlinks to the already-downloaded genomes -------

NCBI_DATA=/shares/sander.imm.uzh/MM/kansasii/data/gtdb_genomes/Mycobacteriaceae/ncbi_dataset/data
LINK_ROOT=/shares/sander.imm.uzh/MM/kansasii/data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii

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
