SPECIES = [
    "kansasii",
    "persicum",
    "pseudokansasii",
    "innocens",
    "attenuatum",
    "ostraviense",
    "gastri",
]

# A fixed, consistent categorical color per species, reused by every figure in
# this repo (main1 and main2) so plots stay comparable across the pipeline.
SPECIES_COLORS = {
    "kansasii": "#1b9e77",
    "persicum": "#d95f02",
    "pseudokansasii": "#7570b3",
    "innocens": "#e7298a",
    "attenuatum": "#66a61e",
    "ostraviense": "#e6ab02",
    "gastri": "#a6761d",
}

REPO_ROOT = __import__("pathlib").Path(__file__).resolve().parent.parent
DATA_ROOT = REPO_ROOT.parent / "data"
GTDB_MLSA_ROOT = DATA_ROOT / "gtdb_genomes" / "Mycobacteriaceae" / "mlsa-kansasii"
SANGER_LINK_DIR = DATA_ROOT / "sanger" / "seq_kansasii"
SCREENING_MAP = DATA_ROOT / "imm" / "screening_map.csv"
