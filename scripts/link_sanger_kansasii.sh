#!/usr/bin/env bash
#
# Create symbolic links to the .ab1 Sanger trace files of all TNR isolates
# listed in data/imm/screening_map.csv, found anywhere under
# data/sanger/seq_mol_routine/***, into a flat directory
# data/sanger/seq_kansasii/.
#
# .fasta/.fa files are intentionally skipped: they can be empty because of
# old base-calling software used historically, so downstream analysis should
# re-basecall from the .ab1 trace instead.
#
# Idempotent: safe to re-run (uses ln -sfn).

set -euo pipefail

ROOT=/shares/sander.imm.uzh/MM/kansasii
SANGER_ROOT="$ROOT/data/sanger/seq_mol_routine"
LINK_DIR="$ROOT/data/sanger/seq_kansasii"
SCREENING_MAP="$ROOT/data/imm/screening_map.csv"

mkdir -p "$LINK_DIR"

TNR_LIST=$(mktemp)
ALL_AB1=$(mktemp)
MATCHED=$(mktemp)
trap 'rm -f "$TNR_LIST" "$ALL_AB1" "$MATCHED"' EXIT

# TNR is column 3; skip header and rows with an empty TNR (reference strains).
awk -F',' 'NR>1 && $3!="" {print $3}' "$SCREENING_MAP" | sort -u > "$TNR_LIST"
n_tnr=$(wc -l < "$TNR_LIST")
echo "Unique TNRs in screening_map.csv: $n_tnr"

find "$SANGER_ROOT" -type f -iname '*.ab1' > "$ALL_AB1"
n_ab1=$(wc -l < "$ALL_AB1")
echo "Total .ab1 files under $SANGER_ROOT: $n_ab1"

# Fixed-string substring match: TNRs are 10-digit numbers, so they cannot
# collide with the 4-digit year/date directory segments in the path.
grep -F -f "$TNR_LIST" "$ALL_AB1" > "$MATCHED" || true
n_matched=$(wc -l < "$MATCHED")
echo "Matched .ab1 files: $n_matched"

count=0
while IFS= read -r filepath; do
    rel="${filepath#"$SANGER_ROOT"/}"
    linkname="${rel//\//__}"
    ln -sfn "$filepath" "$LINK_DIR/$linkname"
    count=$((count + 1))
done < "$MATCHED"

echo "Created/updated $count symlinks in $LINK_DIR"
