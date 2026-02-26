import os
from pathlib import Path

import pandas as pd
import pytest

from nfixplanet.profile import processing


@pytest.mark.parametrize(
    "sample_file,gene_ref,otu_ref",
    [
        (
            os.path.join(
                "tests",
                "references",
                "input",
                "annotate_map",
                "5O_Y1.2_coverage.tsv",
            ),
            os.path.join(
                "tests",
                "references",
                "output",
                "annotate_map",
                "5O_Y1.2_gene_table.tsv",
            ),
            os.path.join(
                "tests",
                "references",
                "output",
                "annotate_map",
                "5O_Y1.2_OTU_table.tsv",
            ),
        ),
    ],
)
def test_annotate_mapping_results(tmp_path, sample_file, gene_ref, otu_ref):
    output_dir = tmp_path / "annotate_map"
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_out = output_dir / Path(gene_ref).name
    otu_out = output_dir / Path(otu_ref).name

    processing.annotate_mapping_results(
        sample_file,
        str(gene_out),
        str(otu_out),
    )

    for produced_path, reference_path in [
        (gene_out, Path(gene_ref)),
        (otu_out, Path(otu_ref)),
    ]:
        assert produced_path.exists()
        assert reference_path.exists()

        produced_df = pd.read_csv(produced_path, sep="\t")
        reference_df = pd.read_csv(reference_path, sep="\t")

        # OTU table requires deterministic sorting
        if "OTU_table" in produced_path.name:
            produced_df = produced_df.sort_values(
                by=list(produced_df.columns)
            ).reset_index(drop=True)

            reference_df = reference_df.sort_values(
                by=list(reference_df.columns)
            ).reset_index(drop=True)
        else:
            produced_df = produced_df.reset_index(drop=True)
            reference_df = reference_df.reset_index(drop=True)

        pd.testing.assert_frame_equal(
            produced_df,
            reference_df,
            check_dtype=False,
        )


@pytest.mark.parametrize(
    "otu_input,ref_dir",
    [
        (
            os.path.join(
                "tests",
                "references",
                "input",
                "profile",
                "OTU_final.tsv",
            ),
            os.path.join(
                "tests",
                "references",
                "output",
                "profile",
            ),
        ),
    ],
)
def test_profile_results(tmp_path, otu_input, ref_dir):
    output_dir = tmp_path / "profile"
    output_dir.mkdir(parents=True, exist_ok=True)

    processing.profile_results(otu_input, str(output_dir))

    ranks = ["d", "p", "c", "o", "f", "g", "s"]

    for rank in ranks:
        produced_path = output_dir / f"OTU_group_summed_by_{rank}.tsv"
        reference_path = Path(ref_dir) / f"OTU_group_summed_by_{rank}.tsv"

        assert produced_path.exists()
        assert reference_path.exists()

        produced_df = (
            pd.read_csv(
                produced_path,
                sep="\t",
                keep_default_na=False,
            )
            .sort_values(by=rank)
            .reset_index(drop=True)
        )

        reference_df = (
            pd.read_csv(
                reference_path,
                sep="\t",
                keep_default_na=False,
            )
            .sort_values(by=rank)
            .reset_index(drop=True)
        )

        pd.testing.assert_frame_equal(
            produced_df,
            reference_df,
            check_dtype=False,
        )
