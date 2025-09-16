import logging
import pandas as pd

logger = logging.getLogger(__name__)


def load_hmm_output(path: str) -> pd.DataFrame:
    column_names = [
        "target_name",
        "accession",
        "query_name",
        "full_accession",
        "full_evalue",
        "full_score",
        "full_bias",
        "domain_evalue",
        "domain_score",
        "domain_bias",
        "exp",
        "reg",
        "clu",
        "ov",
        "env",
        "dom",
        "rep",
        "inc",
    ]
    columns_to_drop = [
        "accession",
        "full_accession",
    ]
    df = pd.read_csv(
        path,
        delim_whitespace=True,
        names=column_names,
        usecols=range(len(column_names)),
        comment="#",
    )
    df = df.drop(columns=columns_to_drop)
    ordered_cols = ["query_name"] + [
        col
        for col in column_names
        if col not in columns_to_drop and col != "query_name"
    ]
    df = df[ordered_cols]
    # These columns are needed for filtering and will be removed before the output
    df[["contig", "gene"]] = df["query_name"].str.rsplit("_", n=1, expand=True)
    return df


def get_best_hmm_hits(df: pd.DataFrame) -> pd.DataFrame:
    # Get the lowest e-value and resolve ties with the highest bit score
    # Resolve ties by taking the first values
    result = df.loc[
        df.groupby("query_name").apply(
            lambda g: g.sort_values(
                ["full_evalue", "full_score"], ascending=[True, False]
            ).index[0]
        )
    ].reset_index(drop=True)

    return result


def extract_nfix(df: pd.DataFrame):
    # Filter by score thresholds
    chl = df[
        ((df["target_name"] == "ChIl") & (df["full_score"] > 194))
        | ((df["target_name"] == "ChlB") & (df["full_score"] > 67))
        | ((df["target_name"] == "ChlN") & (df["full_score"] > 23))
    ]
    
    # Keep only contigs that have all three target_names
    contigs = (
        chl.groupby("contig")["target_name"]
        .apply(lambda x: {"ChIl", "ChlB", "ChlN"}.issubset(set(x)))
        .loc[lambda x: x]
        .index
    )
    logger.debug(list(contigs))

    # Filter ChIl rows within those contigs
    chil = df[(df["target_name"] == "ChIl") & (df["contig"].isin(contigs))]
    chlB = df[(df["target_name"] == "ChlB") & (df["contig"].isin(contigs))]
    chlN = df[(df["target_name"] == "ChlN") & (df["contig"].isin(contigs))]

    logger.debug(chlB)

    return chil



def neighborhood():
    pass


def filter(path: str):
    hmm_output = load_hmm_output(path)
    best_hmm_hits = get_best_hmm_hits(hmm_output)
    extract_nfix(best_hmm_hits)
