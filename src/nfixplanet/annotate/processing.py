import logging
import pandas as pd

from nfixplanet.constants import GeneFamily, GENE_FAMILIES

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
    # NOTE: contig and gene number are needed for neighborhood filtering
    # and will be removed before writing output
    df[["contig", "gene"]] = df["query_name"].str.rsplit("_", n=1, expand=True)
    df["gene"] = df["gene"].astype(int)
    return df


def get_best_hmm_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each query, keep only the best-scoring HMM match:
    lowest e-value wins; ties broken by highest bit score
    """
    return (
        df.sort_values(["full_evalue", "full_score"], ascending=[True, False])
        .groupby("query_name", sort=False)
        .first()
        .reset_index()
    )


def filter_hits_by_score(df: pd.DataFrame, gene_family: GeneFamily) -> pd.DataFrame:
    """
    Keep only rows where the target gene passes its family-specific bitscore threshold.
    """
    mask = pd.Series(False, index=df.index)

    for gene, threshold in gene_family.required.items():
        mask |= (df["target_name"] == gene) & (df["full_score"] > threshold)

    for group in gene_family.alternatives:
        for gene, threshold in group.items():
            mask |= (df["target_name"] == gene) & (df["full_score"] > threshold)

    return df[mask]


def gene_family_to_dict(gf: GeneFamily) -> dict[str, list[int]]:
    result = {}

    result.update({k: [] for k in gf.required})

    if gf.alternatives:
        for alt_dict in gf.alternatives:
            result.update({k: [] for k in alt_dict})

    if gf.optional:
        result.update({k: [] for k in gf.optional})

    return result


def is_within_genomic_context_range(positions: pd.Series, genomic_context_range: int):
    return (positions.max() - positions.min()) <= genomic_context_range


def is_optional_gene_within_genomic_context_range(
    optional_gene_pos: int, positions: pd.Series, genomic_context_range: int
):
    # TODO: ask lucas if optional gene can be min - genomic_context_range
    return positions.max() + genomic_context_range <= optional_gene_pos


def get_valid_contig_indexes(
    hmm_df: pd.DataFrame, gene_family: GeneFamily, genomic_context_range: int
):
    """
    Get indices for rows that meet the following criteria:
    1. all required genes are present
    2. at least one alternative gene is present
    3. all genes in the above set are withint the genomic_context_range
    4. Optional genes are within the genomic_context_range of the largest value
    """
    genes_to_indices = gene_family_to_dict(gene_family)

    for contig, group in hmm_df.groupby("contig"):
        genes_present = set(group["target_name"])
        # NOTE: this will break if there are more than one set of alternatives
        # Fix if that happens
        alt_genes = (
            gene_family.alternatives[0].keys() if gene_family.alternatives else None
        )
        optional_genes = gene_family.optional.keys()

        # All required genes not present
        if not all(g in genes_present for g in gene_family.required):
            continue

        # at least one alternate gene is present
        if alt_genes and not any(g in genes_present for g in alt_genes):
            continue

        # if no duplicates
        if not group["target_name"].duplicated().any():
            if not is_within_genomic_context_range(
                group["gene"], genomic_context_range
            ):
                continue

            postions = group["gene"]

            for idx, row in group.iterrows():
                # combine nifH and vnfH for nif
                if gene_family.name == "nif" and row["target_name"] in alt_genes:
                    genes_to_indices["nifH"].append(idx)
                # combine nifH and vnfH for vnf
                elif gene_family.name == "vnf" and row["target_name"] in alt_genes:
                    genes_to_indices["vnfH"].append(idx)
                # optional genes (nif)
                elif row[
                    "target_name"
                ] in optional_genes and is_optional_gene_within_genomic_context_range(
                    row["gene"], postions, genomic_context_range
                ):
                    genes_to_indices[row["target_name"]].append(idx)
                else:
                    genes_to_indices[row["target_name"]].append(idx)
        else:
            duplicates = group[group["target_name"].duplicated()]
            logger.info(f"Duplicates found for {contig}: {duplicates}")
            # TODO: handle this case later
            continue

    return genes_to_indices


def get_genes_to_df(
    genes_to_indexes: dict[str, list[int]],
    hmm_df: pd.DataFrame,
    gene_family: GeneFamily,
) -> dict[str, pd.DataFrame]:
    genes_to_df = {}
    for gene, indexes in genes_to_indexes.items():
        if not indexes:
            logger.info(f"No hits for {gene} ({gene_family.name})")
            continue
        # logger.debug(f"family: {gene_family.name}\tgene: {gene}")
        # logger.debug(hmm_df.head)
        genes_to_df[gene] = hmm_df.loc[indexes].drop(columns=["contig", "gene"])

    return genes_to_df


def get_filtered_top_hits(
    best_hmm_hits: pd.DataFrame,
    gene_family: GeneFamily,
    genomic_context_range: int,
) -> dict[str, pd.DataFrame]:
    """
    Filter hits that meet the score threshold, gene requirements
    and neighborhood requirements
    """
    # Step 1: filter contigs that are above the bit score threshold
    hits_filtered_by_score = filter_hits_by_score(best_hmm_hits, gene_family)
    if hits_filtered_by_score.empty:
        logger.info(f"No hits above threshold for {gene_family.name}")
        return None

    # Step 2: Get indexes which meet all filtering requirements
    genes_to_indexes = get_valid_contig_indexes(
        hits_filtered_by_score, gene_family, genomic_context_range
    )

    # Step 3: Convert indexes to df
    genes_to_df = get_genes_to_df(genes_to_indexes, best_hmm_hits, gene_family)

    return genes_to_df


def write_tsvs(genes_to_df: dict[str, pd.DataFrame], output_dir: str):
    for gene, df in genes_to_df.items():
        path = f"{output_dir}/{gene}.tsv"
        df.to_csv(path, sep="\t", index=False)
        logger.info(f"Created file: {path}")


def filter_and_write_files(path: str, output_dir: str, genomic_context_range: int):
    hmm_output = load_hmm_output(path)
    best_hmm_hits = get_best_hmm_hits(hmm_output)
    for family in GENE_FAMILIES:
        genes_to_df = get_filtered_top_hits(
            best_hmm_hits, family, genomic_context_range
        )
        if genes_to_df:
            write_tsvs(genes_to_df, output_dir)
