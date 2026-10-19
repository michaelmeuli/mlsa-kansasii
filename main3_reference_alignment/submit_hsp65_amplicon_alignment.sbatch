#!/usr/bin/env bash
#SBATCH --job-name=mlsa-refalign-hsp65
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
ROOT="/shares/sander.imm.uzh/MM/kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

# hsp65 reads against the 72 per-genome hsp65 amplicons (make with scripts/04_make_hsp65_references.py)
python "$REPO_ROOT/main3_reference_alignment/run_reference_alignment.py" \
    --ref-dir "$ROOT/data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons" \
    --out "$ROOT/output/mlsa/main3_reference_alignment_hsp65_amplicons" \
    "$ROOT/data/sanger/seq_kansasii"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/mlsa/main3_reference_alignment_hsp65_amplicons" \
  /shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii submit_hsp65_amplicon_alignment.sbatch /shares/sander.imm.uzh/MM/kansasii/data/sanger/seq_kansasii
