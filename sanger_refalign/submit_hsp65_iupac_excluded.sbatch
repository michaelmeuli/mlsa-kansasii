#!/usr/bin/env bash
#SBATCH --job-name=mlsa-refalign-iupac-excl
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
ROOT="/shares/sander.imm.uzh/MM/kansasii"
OUT="$ROOT/output/mlsa/sanger_refalign_hsp65_iupac"
READS="$ROOT/data/sanger/seq_kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

# F: as C (consensus, every polymorphic site) but without the 3 atypical-hsp65 kansasii genomes
ATYP="GCF_002705785.1 GCF_002705825.1 GCF_002705865.1"
python "$REPO_ROOT/scripts/05_make_iupac_hsp65_references.py" 1 "$OUT/refs_k1_noatypical" $ATYP
python "$REPO_ROOT/sanger_refalign/run_sanger_refalign.py" --iupac \
    --ref-dir "$OUT/refs_k1_noatypical" --out "$OUT/F_consensus_k1_noatypical" "$READS"

bash "$ROOT/repos/kansasii-docs/bin/write_manifest.sh" "$OUT/F_consensus_k1_noatypical" "$REPO_ROOT" submit_hsp65_iupac_excluded.sbatch "$READS"
