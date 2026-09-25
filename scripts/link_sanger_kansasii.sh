#!/usr/bin/env bash
#
# Create symbolic links to the .ab1 Sanger trace files of all TNR isolates
# listed in data/imm/screening_map_link.csv (any of its TNR, TNR_NGS,
# TNR3..TNR6 columns), found anywhere under
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
SCREENING_MAP="$ROOT/data/imm/screening_map_link.csv"

mkdir -p "$LINK_DIR"

TNR_LIST=$(mktemp)
ALL_AB1=$(mktemp)
MATCHED=$(mktemp)
trap 'rm -f "$TNR_LIST" "$ALL_AB1" "$MATCHED"' EXIT

# Locate the TNR columns by header name (their position is not stable across
# versions of the map); skip empty cells (reference strains have no TNR).
awk -F',' '
    NR==1 { for (i=1; i<=NF; i++) if ($i ~ /^TNR(_NGS|[0-9]+)?$/) cols[i]; next }
    { for (i in cols) if ($i != "") print $i }
' "$SCREENING_MAP" | sort -u > "$TNR_LIST"
n_tnr=$(wc -l < "$TNR_LIST")
echo "Unique TNRs in $(basename "$SCREENING_MAP"): $n_tnr"
if grep -qvE '^[0-9]{10}$' "$TNR_LIST"; then
    echo "Non-10-digit TNR in $SCREENING_MAP; refusing to substring-match it" >&2
    exit 1
fi

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
