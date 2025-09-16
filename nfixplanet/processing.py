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


GENE_FAMILIES = {
    # "Chl": {
    #     "ChIl": 194,
    #     "ChlB": 67,
    #     "ChlN": 23,
    # }
    # Core nif
    "nif" : {
        "nifD": 583.6,
        "nifH": 279.5,
        "nifK": 460,
    },
    # nifE and nifN are optional
    "nifE" : {
        "nifD": 583.6,
        "nifH": 279.5,
        "nifK": 460,
        "nifE": 504.6,
    },
    "nifN" : {
        "nifD": 583.6,
        "nifH": 279.5,
        "nifK": 460,
        "nifN": 530,
    },
}

# TODO: add special case for vnf
def write_filtered_file(df: pd.DataFrame, gene_family: str, gene_to_score: dict[str, int], max_dist: int = 10):
    # Only keep the required genes that pass the bit score threshold
    mask = pd.Series(False, index=df.index)
    for gene, threshold in gene_to_score.items():
        mask |= (df["target_name"] == gene) & (df["full_score"] > threshold)
    gene_family_df = df[mask]

    required_genes = set(gene_to_score.keys())

    # Only keep contigs that have all n required genes
    contigs = gene_family_df.groupby("contig")["target_name"].transform(
        lambda x: required_genes.issubset(set(x))
    )
    gene_family_df = gene_family_df[contigs].sort_values(["contig", "gene"])

    # For each contig, find windows of genes that contain all n required genes
    def find_neighborhoods(group):
        # Expand neighborhoods per contig
        results = []
        n = len(group)
        for i in range(n):
            # take window up to max_dist on gene coordinate, not index
            min_gene = group.iloc[i]["gene"]
            window = group[(group["gene"] >= min_gene) & (group["gene"] <= min_gene + max_dist)]
            if required_genes.issubset(set(window["target_name"])):
                results.append(window)
        if results:
            return pd.concat(results)
        return pd.DataFrame(columns=group.columns)

    neighborhood_df = gene_family_df.groupby("contig", group_keys=False).apply(
        find_neighborhoods
    )

    # Subsets for each required gene
    subsets = {gene: neighborhood_df[neighborhood_df["target_name"] == gene]
               for gene in required_genes}
    for gene, subset in subsets.items():
        # Only save nifDKH output without nifEN
        if gene_family in("nifE", "nifN") and gene in ("nifD", "nifH", "nifK"):
            continue
        logger.debug(gene)
        logger.debug(subset)

        # TODO: write files in this function



def filter(path: str):
    hmm_output = load_hmm_output(path)
    best_hmm_hits = get_best_hmm_hits(hmm_output)
    for gene_family, gene_to_score in GENE_FAMILIES.items():
        write_filtered_file(best_hmm_hits, gene_family, gene_to_score)
