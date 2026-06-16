import pytest

from nfixplanet import cli


def test_annotate_command_passes_custom_hmm_format_for_input_hmms(monkeypatch):
    monkeypatch.setattr(cli, "get_version", lambda: "test")
    parser = cli.get_parser()
    args = parser.parse_args(
        [
            "annotate",
            "--input_hmms",
            "input.tbl",
            "--output_directory",
            "out",
            "--custom_hmm_format",
        ]
    )
    called = {}

    def fake_run_annotate(
        input_genomes,
        input_orfs,
        input_hmms,
        output_directory,
        genomic_context_range,
        cpus,
        custom_hmm_format,
    ):
        called["args"] = (
            input_genomes,
            input_orfs,
            input_hmms,
            output_directory,
            genomic_context_range,
            cpus,
            custom_hmm_format,
        )

    monkeypatch.setattr(cli, "run_annotate", fake_run_annotate)

    cli.run_annotate_command(args)

    assert called["args"] == ("", "", "input.tbl", "out", 10, 2, True)


def test_annotate_command_rejects_custom_hmm_format_without_input_hmms(monkeypatch):
    monkeypatch.setattr(cli, "get_version", lambda: "test")
    parser = cli.get_parser()
    args = parser.parse_args(
        [
            "annotate",
            "--input_orfs",
            "input.fna",
            "--output_directory",
            "out",
            "--custom_hmm_format",
        ]
    )

    with pytest.raises(SystemExit, match="--custom_hmm_format can only be used"):
        cli.run_annotate_command(args)


def test_map_command_passes_input_coverage(monkeypatch):
    monkeypatch.setattr(cli, "get_version", lambda: "test")
    parser = cli.get_parser()
    args = parser.parse_args(
        [
            "map",
            "--sample_id",
            "sample-a",
            "--input_coverage",
            "coverm.tsv",
            "--output_directory",
            "out",
        ]
    )
    called = {}

    def fake_run_profile(
        sample_id,
        r1,
        r2,
        single,
        output_directory,
        work_directory,
        cpus,
        input_coverage,
    ):
        called["args"] = (
            sample_id,
            r1,
            r2,
            single,
            output_directory,
            work_directory,
            cpus,
            input_coverage,
        )

    monkeypatch.setattr(cli, "run_profile", fake_run_profile)

    cli.run_profile_command(args)

    assert called["args"] == (
        "sample-a",
        None,
        None,
        None,
        "out",
        "tmp",
        8,
        "coverm.tsv",
    )


def test_map_command_passes_none_for_input_coverage_with_reads(monkeypatch):
    monkeypatch.setattr(cli, "get_version", lambda: "test")
    parser = cli.get_parser()
    args = parser.parse_args(
        [
            "map",
            "--sample_id",
            "sample-a",
            "--read_1",
            "reads_R1.fq.gz",
            "--read_2",
            "reads_R2.fq.gz",
            "--output_directory",
            "out",
        ]
    )
    called = {}

    def fake_run_profile(
        sample_id,
        r1,
        r2,
        single,
        output_directory,
        work_directory,
        cpus,
        input_coverage,
    ):
        called["args"] = (
            sample_id,
            r1,
            r2,
            single,
            output_directory,
            work_directory,
            cpus,
            input_coverage,
        )

    monkeypatch.setattr(cli, "run_profile", fake_run_profile)

    cli.run_profile_command(args)

    assert called["args"] == (
        "sample-a",
        "reads_R1.fq.gz",
        "reads_R2.fq.gz",
        None,
        "out",
        "tmp",
        8,
        None,
    )


def test_map_command_rejects_input_coverage_with_reads(monkeypatch):
    monkeypatch.setattr(cli, "get_version", lambda: "test")
    parser = cli.get_parser()
    args = parser.parse_args(
        [
            "map",
            "--sample_id",
            "sample-a",
            "--input_coverage",
            "coverm.tsv",
            "--single",
            "reads.fq.gz",
            "--output_directory",
            "out",
        ]
    )

    with pytest.raises(SystemExit, match="--input_coverage cannot be used"):
        cli.run_profile_command(args)
