#!/usr/bin/env bash
#SBATCH --job-name=mlsa-refalign-iupac
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
ROOT="/shares/sander.imm.uzh/MM/kansasii"
OUT="$ROOT/output/mlsa/main3_reference_alignment_hsp65_iupac"
AMPL="$ROOT/data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons"
READS="$ROOT/data/sanger/seq_kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

# IUPAC species-consensus references (k1: every polymorphic site, k2: minor allele in >= 2 genomes)
python "$REPO_ROOT/scripts/05_make_iupac_hsp65_references.py" 1 "$OUT/refs_k1"
python "$REPO_ROOT/scripts/05_make_iupac_hsp65_references.py" 2 "$OUT/refs_k2"

RUN="$REPO_ROOT/main3_reference_alignment/run_reference_alignment.py"
# B: IUPAC-aware scoring on the 72 genome amplicons (read-side effect only; A = main3_reference_alignment_hsp65_amplicons)
python "$RUN" --iupac --ref-dir "$AMPL" --out "$OUT/B_amplicons_iupac" "$READS"
# C, D: IUPAC-aware scoring against the species consensus references
python "$RUN" --iupac --ref-dir "$OUT/refs_k1" --out "$OUT/C_consensus_k1" "$READS"
python "$RUN" --iupac --ref-dir "$OUT/refs_k2" --out "$OUT/D_consensus_k2" "$READS"
# E: control, same 7 consensus references with the old (non-IUPAC) scoring
python "$RUN" --ref-dir "$OUT/refs_k2" --out "$OUT/E_consensus_k2_noiupac" "$READS"

bash "$ROOT/repos/kansasii-docs/bin/write_manifest.sh" "$OUT" "$REPO_ROOT" submit_hsp65_iupac_experiment.sbatch "$READS"
