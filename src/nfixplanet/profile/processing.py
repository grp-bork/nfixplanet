import re
import pandas as pd
from pathlib import Path

from nfixplanet.constants import (
    MAPPING_REFERENCE,
    TARGET_ANNOTATIONS,
    TAXONOMY_REFERENCE,
)


def load_otu_table(sample_path: str) -> pd.DataFrame:
    """Get means for contigs"""
    df = pd.read_csv(sample_path, sep="\t")

    # If paired and single
    if df.shape[1] == 5:
        mean_cols = df.columns[df.columns.str.contains("Mean")]
        row_sum = df[mean_cols].sum(axis=1)
        df["Mean"] = row_sum
    # Just rename columns
    else:
        df.columns = ["Contig", "Covered_Bases", "Mean"]

    # only return mean
    return df[["Contig", "Mean"]]


def extract_sample_name(path: str) -> str:
    match = re.search(r"(.*?)_coverage\.tsv", Path(path).name)
    if not match:
        raise ValueError("Could not extract sample name")
    return match.group(1)


def compute_target_otu_table(merged_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    """Filter merged data to target annotations and rename Mean column to sample name."""
    subset = merged_df[merged_df["annotation"].isin(TARGET_ANNOTATIONS)]
    otu_table = subset[["map_id", "annotation", "Mean"]].copy()
    otu_table.rename(columns={"Mean": sample_name}, inplace=True)
    return otu_table


def compute_genome_normalized_table(
    reference_df: pd.DataFrame, otu_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Compute genome-level normalized abundances by averaging values across genes per genome_contig.
    """

    # attach genome_contig
    contig_map = reference_df[["map_id", "genome_contig"]]
    merged = contig_map.merge(otu_df, on="map_id", how="inner")

    # sum per genome_contig
    grouped = (
        merged.drop(columns=["map_id", "annotation"])
        .groupby("genome_contig", as_index=False)
        .sum()
    )

    # count genes per genome
    gene_counts = merged.groupby("genome_contig").size().reset_index(name="Freq")

    normalized = grouped.merge(gene_counts, on="genome_contig")

    sample_cols = normalized.columns.difference(["genome_contig", "Freq"])

    normalized[sample_cols] = normalized[sample_cols].div(normalized["Freq"], axis=0)

    normalized = normalized.drop(columns=["Freq"])

    # Sort output to match R
    normalized = normalized.sort_values(
        "genome_contig", key=lambda x: x.str.lower()
    ).reset_index(drop=True)

    return normalized


def annotate_mapping_results(coverage_file: str, gene_file: str, otu_file: str):
    """Process coverage data to generate gene-level and genome-normalized OTU tables."""
    reference_df = pd.read_csv(MAPPING_REFERENCE, sep="\t")
    otu_df = load_otu_table(coverage_file)

    sample_name = extract_sample_name(coverage_file)

    # merge with annotation
    merged_df = reference_df.merge(
        otu_df, left_on="map_id", right_on="Contig", how="inner"
    )

    # Aggregate Mean coverage values by annotation and label column with sample name
    gene_result = (
        merged_df.groupby("annotation", as_index=False)["Mean"]
        .sum()
        .rename(columns={"Mean": sample_name})
    )

    otu_table = compute_target_otu_table(merged_df, sample_name)
    otu_normalized = compute_genome_normalized_table(reference_df, otu_table)

    gene_result.to_csv(gene_file, sep="\t", index=False, float_format="%.15g")

    otu_normalized.to_csv(otu_file, sep="\t", index=False, float_format="%.15g")


def sum_otu_by_group(
    group_col: str,
    otu_df: pd.DataFrame,
    tax_df: pd.DataFrame,
    output_dir: str,
    out_prefix: str = "OTU_group_summed",
):
    """Sum OTU abundance values grouped by a specified taxonomy column and write to file."""
    tax_map = tax_df[["genome_contig", group_col]]
    merged = otu_df.merge(tax_map, on="genome_contig", how="inner")

    numeric_cols = merged.select_dtypes(include="number").columns.tolist()

    grouped = merged.groupby(group_col, dropna=False)[numeric_cols].sum().reset_index()

    out_file = f"{output_dir}/{out_prefix}_by_{group_col}.tsv"
    grouped.to_csv(out_file, sep="\t", index=False, float_format="%.15g")


def profile_results(otu_path: str, output_dir: str):
    """Generate taxonomy-level OTU summaries for standard taxonomic ranks."""
    # Stop autoconverting "" to "NaN"
    tax_df = pd.read_csv(TAXONOMY_REFERENCE, sep="\t", keep_default_na=False)
    otu_df = pd.read_csv(otu_path, sep="\t")

    for rank in ["d", "p", "c", "o", "f", "g", "s"]:
        sum_otu_by_group(rank, otu_df, tax_df, output_dir)
