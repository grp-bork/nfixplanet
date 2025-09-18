import logging
import pandas as pd

from .config import GeneFamily, GENE_FAMILIES

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
    df = pd.read_csv(
        path,
        delim_whitespace=True,
        names=column_names,
        usecols=range(len(column_names)),
        comment="#",
    )
    ordered_cols = ["query_name", "accession"] + [
        col for col in column_names if col not in ["query_name", "accession"]
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


def filter_top_hits_by_genes(
    df: pd.DataFrame, gene_family: GeneFamily, max_dist: int = 10
) -> dict[str, pd.DataFrame] | None:
    required_genes = gene_family.required
    alternative_groups = gene_family.alternatives

    # Build filter: any gene in required and at least one in alternatives above threshold
    mask = pd.Series(False, index=df.index)

    for gene, threshold in required_genes.items():
        mask |= (df["target_name"] == gene) & (df["full_score"] > threshold)

    if alternative_groups:
        for group in alternative_groups:
            for gene, threshold in group.items():
                mask |= (df["target_name"] == gene) & (df["full_score"] > threshold)

    gene_family_df = df[mask]

    if gene_family_df.empty:
        logger.info(f"No hits found for {gene_family.name}")
        return None

    def contains_required_gene_combination(genes: set[str]) -> bool:
        # all required must be present
        if not required_genes.keys() <= genes:
            return False
        # each alt-group: at least one must be present
        if alternative_groups:
            for group in alternative_groups:
                if genes.isdisjoint(group.keys()):
                    return False
        return True

    # Only keep contigs that have all n required genes
    contigs = gene_family_df.groupby("contig")["target_name"].transform(
        lambda x: contains_required_gene_combination(set(x))
    )
    logger.debug(gene_family_df)
    gene_family_df = gene_family_df[contigs].sort_values(["contig", "gene"])

    # For each contig, find windows of genes that contain all n required genes
    def find_neighborhoods(group):
        # Expand neighborhoods per contig
        results = []
        n = len(group)
        for i in range(n):
            # take window up to max_dist on gene coordinate, not index
            min_gene = group.iloc[i]["gene"]
            window = group[
                (group["gene"] >= min_gene) & (group["gene"] <= min_gene + max_dist)
            ]
            if contains_required_gene_combination(set(window["target_name"])):
                results.append(window)
        if results:
            return pd.concat(results)
        return pd.DataFrame(columns=group.columns)

    neighborhood_df = gene_family_df.groupby("contig", group_keys=False).apply(
        find_neighborhoods
    )

    # Subsets for each required gene
    subsets: dict[str, pd.DataFrame] = {}
    for gene in required_genes:
        subsets[gene] = neighborhood_df[neighborhood_df["target_name"] == gene]

    if alternative_groups:
        for group in alternative_groups:
            for gene in group:
                subsets[gene] = neighborhood_df[neighborhood_df["target_name"] == gene]
    return subsets


def write_tsv(
    genes_to_hits: dict[str, pd.DataFrame], gene_family_name: str, output_dir: str
):
    for gene, subset in genes_to_hits.items():
        # Only save nifDKH output without nifEN
        if gene_family_name in ("nifE", "nifN") and gene in ("nifD", "nifH", "nifK"):
            continue
        if subset.empty:
            continue
        # logger.debug(f"family: {gene_family_name}\tgene: {gene}")
        # logger.debug(subset.head)
        path = f"{output_dir}/{gene}.tsv"
        subset = subset.drop(["contig", "gene"], axis=1)
        subset.to_csv(path, sep="\t", index=False)


def filter_and_write_files(path: str, output_dir: str):
    hmm_output = load_hmm_output(path)
    best_hmm_hits = get_best_hmm_hits(hmm_output)
    for family in GENE_FAMILIES:
        genes_to_hits = filter_top_hits_by_genes(best_hmm_hits, family)
        if genes_to_hits:
            write_tsv(genes_to_hits, family.name, output_dir)
