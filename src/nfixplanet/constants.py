from pathlib import Path
from dataclasses import dataclass
from importlib import resources

DATA_DIR = resources.files("nfixplanet.reference_data")

######################################
# Required by Annotate pipeline
######################################
HMM_PROFILE_PATH = str(DATA_DIR / "hmm_profiles" / "nfixplanet_models.hmm")


@dataclass
class GeneFamily:
    name: str
    required: dict[str, float]
    alternatives: list[dict[str, float]] | None
    optional: dict[str, float] | None


GENE_FAMILIES = [
    GeneFamily(
        name="nif",
        required={"nifD": 583.6, "nifK": 460.0},
        alternatives=[{"nifH": 279.5, "vnfH": 150.6}],
        optional = {"nifE": 504.6, "nifN": 530},
    ),
    GeneFamily(
        name="vnf",
        required={"vnfD": 933.0, "vnfK": 1151.5},
        alternatives=[{"nifH": 279.5, "vnfH": 150.6}],
        optional={},
    ),
    GeneFamily(
        name="anf",
        required={"anfD": 866.0, "anfK": 1016.1, "anfH": 626.6},
        alternatives=[],
        optional={},
    ),
    GeneFamily(
        name="chl",
        required={"ChlL": 194.0, "ChlB": 67.0, "ChlN": 23.0},
        alternatives=[],
        optional={},
    ),
    GeneFamily(
        name="nfl",
        required={"nflD": 88.0, "nflH": 509.0},
        alternatives=[],
        optional={},
    ),
]


ANNOTATE_TOOLS = [
    "prodigal",
    "hmmscan",
    "zcat",
]

######################################
# Required by Map pipeline
######################################

MAP_TOOLS = [
    "wget",
    "gunzip",
    "fastp",
    "hostile",
    "coverm",
]

# TODO: Make cache dir changable
HOSTILE_CACHE_DIR = Path("~/.local/share/hostile").expanduser()

NFIXPLANET_CACHE_DIR = Path("~/.local/share/nfixplanet").expanduser()
REFERENCE_FASTA_URL = "https://zenodo.org/records/18772121/files/nfix_reference_v1.0.0.fna.gz?download=1"
REFERENCE_FASTA_NAME = NFIXPLANET_CACHE_DIR / "nfix_reference_v1.0.0.fna"
REFERENCE_INDEX_NAME = NFIXPLANET_CACHE_DIR / "nfix_reference_v1.0.0.fna.mmi"

MAPPING_REFERENCE = str(DATA_DIR / "map" / "mapping_ref_22.01.26.tsv")
TAXONOMY_REFERENCE = str(DATA_DIR / "map" / "taxonomy_ref_27.01.26.tsv")


TARGET_ANNOTATIONS = {
    "nifD",
    "nifK",
    "anfD",
    "anfK",
    "vnfD",
    "vnfK",
}
