import re
import pandas as pd
from pathlib import Path

from nfixplanet.constants import (
    MAPPING_REFERENCE,
    TARGET_ANNOTATIONS,
    TAXONOMY_REFERENCE,
)


def load_reference(reference_path: str) -> pd.DataFrame:
    return pd.read_csv(reference_path, sep="\t")


def load_otu_table(sample_path: str) -> pd.DataFrame:
    df = pd.read_csv(sample_path, sep="\t")

    # paired reads
    if df.shape[1] == 5:
        df["total_mean"] = (
            df["nfix_reference_final.mmi.R1.clean_1.fastq.gz.Mean"]
            + df["nfix_reference_final.mmi.RS.clean.fastq.gz.Mean"]
        )
        df = df[["Contig", "total_mean"]]
        df.columns = ["Contig", "Mean"]
    else:
        df.columns = ["Contig", "Covered_Bases", "Mean"]
        df = df[["Contig", "Mean"]]

    return df


def extract_sample_name(path: str) -> str:
    match = re.search(r"(.*?)_sample_coverage\.tsv", Path(path).name)
    if not match:
        raise ValueError("Could not extract sample name")
    return match.group(1)


def merge_with_annotation(
    reference_df: pd.DataFrame, otu_df: pd.DataFrame
) -> pd.DataFrame:
    return reference_df.merge(otu_df, left_on="map_id", right_on="Contig", how="inner")


def compute_gene_level(merged_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    gene_result = (
        merged_df.groupby("annotation", as_index=False)["Mean"]
        .sum()
        .rename(columns={"Mean": sample_name})
    )
    return gene_result


def compute_target_otu_table(merged_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    subset = merged_df[merged_df["annotation"].isin(TARGET_ANNOTATIONS)]
    otu_result = subset[["map_id", "annotation", "Mean"]].copy()
    otu_result.rename(columns={"Mean": sample_name}, inplace=True)
    return otu_result


def compute_genome_normalized_table(
    reference_df: pd.DataFrame, otu_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Equivalent of:
    - merge genome_contig
    - sum per genome
    - divide by gene count per genome
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


def annotate_mapping_results(sample_file: str, output_dir: str):

    reference_df = load_reference(MAPPING_REFERENCE)
    otu_df = load_otu_table(sample_file)

    sample_name = extract_sample_name(sample_file)

    merged_df = merge_with_annotation(reference_df, otu_df)

    gene_result = compute_gene_level(merged_df, sample_name)
    otu_result = compute_target_otu_table(merged_df, sample_name)
    otu_normalized = compute_genome_normalized_table(reference_df, otu_result)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    gene_result.to_csv(
        output_path / "gene_final.tsv", sep="\t", index=False, float_format="%.15g"
    )

    otu_result.to_csv(
        output_path / "nifDK_final.tsv", sep="\t", index=False, float_format="%.15g"
    )

    otu_normalized.to_csv(
        output_path / "OTU_final.tsv", sep="\t", index=False, float_format="%.15g"
    )


def sum_otu_by_group(
    group_col: str,
    otu_df: pd.DataFrame,
    tax_df: pd.DataFrame,
    output_dir: str,
    out_prefix: str = "OTU_group_summed",
):
    tax_map = tax_df[["genome_contig", group_col]]
    merged = otu_df.merge(tax_map, on="genome_contig", how="inner")

    numeric_cols = merged.select_dtypes(include="number").columns.tolist()

    grouped = merged.groupby(group_col, dropna=False)[numeric_cols].sum().reset_index()

    out_file = f"{output_dir}/{out_prefix}_by_{group_col}.tsv"
    grouped.to_csv(out_file, sep="\t", index=False, float_format="%.15g")


def profile_results(otu_path: str, output_dir: str):
    # Stop autoconverting "" to "NaN"
    tax_df = pd.read_csv(TAXONOMY_REFERENCE, sep="\t", keep_default_na=False)
    otu_df = pd.read_csv(otu_path, sep="\t")

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for rank in ["d", "p", "c", "o", "f", "g", "s"]:
        sum_otu_by_group(rank, otu_df, tax_df, output_dir)
