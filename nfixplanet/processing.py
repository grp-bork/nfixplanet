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
    chl = df[
        ((df["target_name"] == "ChIl") & (df["full_score"] > 194))
        | ((df["target_name"] == "ChlB") & (df["full_score"] > 67))
        | ((df["target_name"] == "ChlN") & (df["full_score"] > 23))
    ]

    # Keep only contigs that have all three target_names
    contigs = (
        chl.groupby("contig")["target_name"]
        .apply(lambda x: {"ChIl", "ChlB", "ChlN"}.issubset(set(x)))
    )
    valid_contigs = contigs[contigs].index

    # Filter to only those contigs
    chl = chl[chl["contig"].isin(valid_contigs)].sort_values("gene")

    targets = ["ChlN", "ChlB", "ChIl"]
    targets_to_gene = {gene: None for gene in targets}
    
    neighborhood_queries = set()
    for index, row in chl.iterrows():
        # valid contig
        if row["query_name"] in neighborhood_queries:
            continue
        targets_to_gene[row["target_name"]] = row["gene"]
        # all target genes not found yet
        if not all(targets_to_gene.values()):
            continue
        max_gene = max(targets_to_gene.values())
        min_gene = min(targets_to_gene.values())

        if max_gene - min_gene <= max_dist:
            for gene in targets_to_gene.values():
                # build the correct query_name
                neighborhood_queries.add(row["contig"] + "_" + str(gene))

    neighborhood_df = chl[chl["query_name"].isin(neighborhood_queries)]


    chIl = neighborhood_df[neighborhood_df["target_name"] == "ChIl"]
    chlB = neighborhood_df[neighborhood_df["target_name"] == "ChlB"]
    chlN = neighborhood_df[neighborhood_df["target_name"] == "ChlN"]
    logger.debug(chIl)
    logger.debug(chlB)
    logger.debug(chlN)

    # Final filter
    return None



def neighborhood():
    pass


def filter(path: str):
    hmm_output = load_hmm_output(path)
    best_hmm_hits = get_best_hmm_hits(hmm_output)
    extract_nfix(best_hmm_hits)
    #logger.debug(x)
