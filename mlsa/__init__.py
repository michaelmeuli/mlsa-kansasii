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
DATA_ROOT = REPO_ROOT.parent.parent / "data"
GTDB_MLSA_ROOT = DATA_ROOT / "gtdb_genomes" / "Mycobacteriaceae" / "mlsa-kansasii"
SANGER_LINK_DIR = DATA_ROOT / "sanger" / "seq_kansasii"
SCREENING_MAP = DATA_ROOT / "imm" / "screening_map_link.csv"
# Every column of SCREENING_MAP holding a TNR of the row's isolate: TNR is
# the primary one, the others are further samples of the same isolate.
SCREENING_MAP_TNR_COLUMNS = ("TNR", "TNR_NGS", "TNR3", "TNR4", "TNR5", "TNR6")
