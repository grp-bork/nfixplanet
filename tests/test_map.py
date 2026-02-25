import os
import pandas as pd
from pathlib import Path

from nfixplanet.map import processing


def test_annotate_mapping_results(tmp_path):
    # input sample file
    sample_file = os.path.join(
        "tests",
        "references",
        "input",
        "annotate_map",
        "5O_Y1.2_sample_coverage.tsv",
    )

    # run function with temporary output dir
    output_dir = tmp_path / "annotate_map"
    processing.annotate_mapping_results(sample_file, str(output_dir))

    # expected output files
    expected_files = {
        "gene_final.tsv": "gene_final.tsv",
        "nifDK_final.tsv": "nifDK_final.tsv",
        "OTU_final.tsv": "OTU_final_sorted.tsv",  # reference is sorted
    }

    for produced_name, reference_name in expected_files.items():
        produced_path = output_dir / produced_name
        reference_path = Path(
            "tests",
            "references",
            "output",
            "annotate_map",
            reference_name,
        )

        assert produced_path.exists(), f"Missing output file: {produced_name}"
        assert reference_path.exists(), f"Missing reference file: {reference_name}"

        produced_df = pd.read_csv(produced_path, sep="\t")
        reference_df = pd.read_csv(reference_path, sep="\t")

        # For OTU table, ensure deterministic ordering before compare
        if produced_name == "OTU_final.tsv":
            produced_df = produced_df.sort_values(
                by=list(produced_df.columns)
            ).reset_index(drop=True)
            reference_df = reference_df.sort_values(
                by=list(reference_df.columns)
            ).reset_index(drop=True)
        else:
            produced_df = produced_df.reset_index(drop=True)
            reference_df = reference_df.reset_index(drop=True)

        produced_values = produced_df.values
        reference_values = reference_df.values

        assert (
            produced_values == reference_values
        ).all(), f"Mismatch in {produced_name}"


def test_profile_results(tmp_path):
    # input OTU table
    otu_path = os.path.join(
        "tests",
        "references",
        "input",
        "profile",
        "OTU_final.tsv",
    )

    output_dir = tmp_path / "profile"
    processing.profile_results(otu_path, str(output_dir))

    ranks = ["d", "p", "c", "o", "f", "g", "s"]

    for rank in ranks:
        produced_path = output_dir / f"OTU_group_summed_by_{rank}.tsv"
        reference_path = Path(
            "tests",
            "references",
            "output",
            "profile",
            f"OTU_group_summed_by_{rank}.tsv",
        )

        assert produced_path.exists(), f"Missing output for rank {rank}"
        assert reference_path.exists(), f"Missing reference for rank {rank}"


        produced_df = pd.read_csv(produced_path, sep="\t", keep_default_na=False ).reset_index(drop=True)
        reference_df = pd.read_csv(reference_path, sep="\t", keep_default_na=False).reset_index(drop=True)
        

        # sort rows by rank
        produced_df = produced_df.sort_values(by=rank).reset_index(drop=True)
        reference_df = reference_df.sort_values(by=rank).reset_index(drop=True)

        pd.testing.assert_frame_equal(
            produced_df,
            reference_df,
            check_dtype=False,
        )