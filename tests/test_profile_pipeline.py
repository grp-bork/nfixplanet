from pathlib import Path

import pandas as pd

from nfixplanet.profile.pipeline import run_profile


def test_profile_pipeline_with_input_coverage(tmp_path):
    sample_id = "5O_Y1.2"
    output_dir = tmp_path / "output_annotate_map"
    input_coverage = Path(
        "tests/references/input/annotate_map/5O_Y1.2_coverage.tsv"
    )
    reference_dir = Path("tests/references/output/annotate_map")
    reference_profile_dir = Path("tests/references/output/profile")

    run_profile(
        sample_id,
        None,
        None,
        None,
        str(output_dir),
        str(tmp_path / "work"),
        8,
        str(input_coverage),
    )

    for filename in [
        "5O_Y1.2_gene_table.tsv",
        "5O_Y1.2_OTU_table.tsv",
    ]:
        produced_path = output_dir / filename
        reference_path = reference_dir / filename

        assert produced_path.exists()

        produced_df = pd.read_csv(produced_path, sep="\t")
        reference_df = pd.read_csv(reference_path, sep="\t")

        if "OTU_table" in filename:
            produced_df = produced_df.sort_values(
                by=list(produced_df.columns)
            ).reset_index(drop=True)
            reference_df = reference_df.sort_values(
                by=list(reference_df.columns)
            ).reset_index(drop=True)

        pd.testing.assert_frame_equal(
            produced_df,
            reference_df,
            check_dtype=False,
        )

    for rank in ["d", "p", "c", "o", "f", "g", "s"]:
        filename = f"OTU_group_summed_by_{rank}.tsv"
        produced_path = output_dir / filename
        reference_path = reference_profile_dir / filename

        assert produced_path.exists()

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

    assert not (output_dir / "5O_Y1.2_coverage.tsv").exists()
