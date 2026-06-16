from pathlib import Path

import pandas as pd

from nfixplanet.profile import pipeline
from nfixplanet.profile.pipeline import run_profile


def test_profile_pipeline_with_reads_orchestrates_mapping(monkeypatch, tmp_path):
    calls = []
    output_dir = tmp_path / "output"
    work_dir = tmp_path / "work"
    hostile_cache_dir = tmp_path / "hostile_cache"
    minimap_index = tmp_path / "reference.mmi"

    hostile_cache_dir.mkdir()
    (hostile_cache_dir / "existing-index").write_text("cached")
    minimap_index.write_text("indexed")

    monkeypatch.setattr(pipeline, "HOSTILE_CACHE_DIR", hostile_cache_dir)
    monkeypatch.setattr(pipeline, "REFERENCE_INDEX_NAME", minimap_index)
    monkeypatch.setattr(
        pipeline.utils,
        "check_external_tools",
        lambda tools: calls.append(("check_external_tools", tools)),
    )
    monkeypatch.setattr(
        pipeline.tools,
        "check_files_exist",
        lambda paths: calls.append(("check_files_exist", paths)),
    )
    monkeypatch.setattr(
        pipeline.tools,
        "fastp_paired",
        lambda **kwargs: calls.append(("fastp_paired", kwargs)),
    )
    monkeypatch.setattr(
        pipeline.tools,
        "hostile_index_fetch",
        lambda: calls.append(("hostile_index_fetch",)),
    )
    monkeypatch.setattr(
        pipeline.tools,
        "hostile_clean_paired",
        lambda r1, r2, out_dir, cpus: (
            calls.append(
                (
                    "hostile_clean_paired",
                    {
                        "r1": r1,
                        "r2": r2,
                        "out_dir": out_dir,
                        "cpus": cpus,
                    },
                )
            )
            or (f"{out_dir}/R1.clean_1.fastq.gz", f"{out_dir}/R2.clean_2.fastq.gz")
        ),
    )
    monkeypatch.setattr(
        pipeline.tools,
        "build_minimap_index",
        lambda cpus: calls.append(("build_minimap_index", cpus)),
    )
    monkeypatch.setattr(
        pipeline.tools,
        "coverm_contig",
        lambda **kwargs: calls.append(("coverm_contig", kwargs)),
    )
    monkeypatch.setattr(
        pipeline,
        "process_coverage_results",
        lambda sample_id, coverage_file, output_dir: calls.append(
            (
                "process_coverage_results",
                {
                    "sample_id": sample_id,
                    "coverage_file": coverage_file,
                    "output_dir": output_dir,
                },
            )
        ),
    )

    run_profile(
        "sample-a",
        "reads_R1.fq.gz",
        "reads_R2.fq.gz",
        None,
        str(output_dir),
        str(work_dir),
        8,
        None,
    )

    assert output_dir.exists()
    assert (work_dir / "processed_fastqs").exists()
    assert (work_dir / "cleaned_fastqs").exists()
    assert ("hostile_index_fetch",) not in calls
    assert ("build_minimap_index", 3) not in calls

    assert calls == [
        ("check_external_tools", pipeline.MAP_TOOLS),
        ("check_files_exist", ["reads_R1.fq.gz", "reads_R2.fq.gz"]),
        (
            "fastp_paired",
            {
                "r1_in": "reads_R1.fq.gz",
                "r2_in": "reads_R2.fq.gz",
                "r1_out": f"{work_dir}/processed_fastqs/R1.fq.gz",
                "r2_out": f"{work_dir}/processed_fastqs/R2.fq.gz",
                "unpaired": f"{work_dir}/processed_fastqs/unpaired.fq",
                "html": f"{work_dir}/processed_fastqs/paired.html",
                "json": f"{work_dir}/processed_fastqs/paired.json",
                "cpus": 4,
            },
        ),
        (
            "hostile_clean_paired",
            {
                "r1": f"{work_dir}/processed_fastqs/R1.fq.gz",
                "r2": f"{work_dir}/processed_fastqs/R2.fq.gz",
                "out_dir": f"{work_dir}/cleaned_fastqs",
                "cpus": 8,
            },
        ),
        (
            "coverm_contig",
            {
                "r1": f"{work_dir}/cleaned_fastqs/R1.clean_1.fastq.gz",
                "r2": f"{work_dir}/cleaned_fastqs/R2.clean_2.fastq.gz",
                "single": None,
                "reference_index": str(minimap_index),
                "output_file": f"{output_dir}/sample-a_coverage.tsv",
                "cpus": 8,
            },
        ),
        (
            "process_coverage_results",
            {
                "sample_id": "sample-a",
                "coverage_file": f"{output_dir}/sample-a_coverage.tsv",
                "output_dir": str(output_dir),
            },
        ),
    ]


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
