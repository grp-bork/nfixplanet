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
