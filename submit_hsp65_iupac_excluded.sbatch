#!/usr/bin/env bash
# IUPAC experiment F: consensus without the atypical-hsp65 genomes -> output/mlsa-kansasii/results/sanger_refalign_hsp65_iupac/F_consensus_k1_noatypical/
#SBATCH --job-name=mlsa-refalign-iupac-excl
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa-kansasii/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUT_MLSA/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

OUT="$OUTPUT/sanger_refalign_hsp65_iupac"
# F: as C (consensus, every polymorphic site) but without the 3 atypical-hsp65 kansasii genomes (ATYPICAL_HSP65_GENOMES)
ATYP="GCF_002705785.1 GCF_002705825.1 GCF_002705865.1"
python "$REPO/scripts/07_make_iupac_hsp65_references.py" 1 "$OUT/refs_k1_noatypical" $ATYP
python "$REPO/scripts/09_sanger_refalign.py" --iupac --ref-dir "$OUT/refs_k1_noatypical" \
    --out "$OUT/F_consensus_k1_noatypical" "$SANGER_LINK_DIR"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUT/F_consensus_k1_noatypical" "$REPO" submit_hsp65_iupac_excluded.sbatch "$SANGER_LINK_DIR"
