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
    df["gene"] = df["gene"].astype(int)
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


def extract_nfix(df: pd.DataFrame, max_dist: int = 10):
    # Filter by score thresholds
    gene_family = df[
        ((df["target_name"] == "ChIl") & (df["full_score"] > 194))
        | ((df["target_name"] == "ChlB") & (df["full_score"] > 67))
        | ((df["target_name"] == "ChlN") & (df["full_score"] > 23))
    ]

    # Keep only contigs that have all three target_names
    contigs = (
        gene_family.groupby("contig")["target_name"]
        .transform(lambda x: {"ChIl", "ChlB", "ChlN"}.issubset(set(x)))
    )
    gene_family = gene_family[contigs].sort_values(["contig", "gene"])

    # For each contig, find windows of genes that contain all three targets
    def find_neighborhoods(g):
        # Expand neighborhoods per contig
        results = []
        n = len(g)
        for i in range(n):
            # take window up to max_dist on gene coordinate, not index
            min_gene = g.iloc[i]["gene"]
            window = g[(g["gene"] >= min_gene) & (g["gene"] <= min_gene + max_dist)]
            if {"ChIl", "ChlB", "ChlN"}.issubset(set(window["target_name"])):
                results.append(window)
        if results:
            return pd.concat(results)
        return pd.DataFrame(columns=g.columns)

    neighborhood_df = gene_family.groupby("contig", group_keys=False).apply(find_neighborhoods)

    # Subsets
    chIl = neighborhood_df[neighborhood_df["target_name"] == "ChIl"]
    chlB = neighborhood_df[neighborhood_df["target_name"] == "ChlB"]
    chlN = neighborhood_df[neighborhood_df["target_name"] == "ChlN"]
    logger.debug(chIl)
    logger.debug(chlB)
    logger.debug(chlN)

    return neighborhood_df



def filter(path: str):
    hmm_output = load_hmm_output(path)
    best_hmm_hits = get_best_hmm_hits(hmm_output)
    extract_nfix(best_hmm_hits)
    #logger.debug(x)
