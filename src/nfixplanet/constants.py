from dataclasses import dataclass
from importlib import resources

DATA_DIR = resources.files("nfixplanet.reference_data.hmm_profiles")
HMM_PROFILE_PATH = str(DATA_DIR / "nfixplanet_models.hmm")

@dataclass
class GeneFamily:
    name: str
    required: dict[str, float]
    alternatives: list[dict[str, float]]


GENE_FAMILIES = [
    GeneFamily(
        name="nif",
        required={"nifD": 583.6, "nifK": 460.0},
        alternatives=[{"nifH": 279.5, "vnfH": 150.6}],
    ),
    # NOTE: nifN and nifE are optional. It was easier to
    # add them as a different gene fmaily instead of allowing
    # optional outputs
    GeneFamily(
        name="nifE",
        required={"nifD": 583.6, "nifK": 460.0, "nifE": 504.6},
        alternatives=[{"nifH": 279.5, "vnfH": 150.6}],
    ),
    GeneFamily(
        name="nifN",
        required={"nifD": 583.6, "nifK": 460.0, "nifN": 530},
        alternatives=[{"nifH": 279.5, "vnfH": 150.6}],
    ),
    GeneFamily(
        name="vnf",
        required={"vnfD": 933.0, "vnfK": 1151.5},
        alternatives=[{"nifH": 279.5, "vnfH": 150.6}],
    ),
    GeneFamily(
        name="anf",
        required={"anfD": 866.0, "anfK": 1016.1, "anfH": 626.6},
        alternatives=[],
    ),
    GeneFamily(
        name="chl",
        required={"ChIl": 194.0, "ChlB": 67.0, "ChlN": 23.0},
        alternatives=[],
    ),
    GeneFamily(
        name="nfl",
        required={"nflD": 88.0, "nflH": 509.0},
        alternatives=[],
    ),
]
