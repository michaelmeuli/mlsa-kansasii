#!/usr/bin/env bash
# IUPAC experiment A-E (README) -> output/mlsa-kansasii/results/sanger_refalign_hsp65_iupac/
#SBATCH --job-name=mlsa-refalign-iupac
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa-kansasii/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUT_MLSA/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

OUT="$OUTPUT/sanger_refalign_hsp65_iupac"
RUN="$REPO/scripts/09_sanger_refalign.py"

# IUPAC species-consensus references (k1: every polymorphic site, k2: minor allele in >= 2 genomes)
python "$REPO/scripts/07_make_iupac_hsp65_references.py" 1 "$OUT/refs_k1"
python "$REPO/scripts/07_make_iupac_hsp65_references.py" 2 "$OUT/refs_k2"

# B: IUPAC-aware scoring on the 72 genome amplicons (read-side effect only; A = sanger_refalign_hsp65_amplicons)
python "$RUN" --iupac --ref-dir "$HSP65_AMPLICONS" --out "$OUT/B_amplicons_iupac" "$SANGER_LINK_DIR"
# C, D: IUPAC-aware scoring against the species consensus references
python "$RUN" --iupac --ref-dir "$OUT/refs_k1" --out "$OUT/C_consensus_k1" "$SANGER_LINK_DIR"
python "$RUN" --iupac --ref-dir "$OUT/refs_k2" --out "$OUT/D_consensus_k2" "$SANGER_LINK_DIR"
# E: control, same 7 consensus references with the old (non-IUPAC) scoring
python "$RUN" --ref-dir "$OUT/refs_k2" --out "$OUT/E_consensus_k2_noiupac" "$SANGER_LINK_DIR"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/sanger_refalign_hsp65_iupac" "$REPO" submit_hsp65_iupac_experiment.sbatch "$SANGER_LINK_DIR" "$HSP65_AMPLICONS"
