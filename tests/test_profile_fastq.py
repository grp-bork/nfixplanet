import gzip
from pathlib import Path

from nfixplanet.profile import fastq


def test_prepare_fastqs_promotes_unpaired_reads_to_single_reads(monkeypatch, tmp_path):
    r1 = tmp_path / "reads_R1.fq"
    r2 = tmp_path / "reads_R2.fq"
    fastp_dir = tmp_path / "fastp"
    unpaired_content = b"@unpaired\nACGT\n+\n!!!!\n"
    calls = []

    r1.write_text("@r1\nACGT\n+\n!!!!\n")
    r2.write_text("@r2\nTGCA\n+\n!!!!\n")
    fastp_dir.mkdir()

    def fake_fastp_paired(**kwargs):
        calls.append(kwargs)
        Path(kwargs["unpaired"]).write_bytes(unpaired_content)

    monkeypatch.setattr(fastq.tools, "fastp_paired", fake_fastp_paired)

    processed_reads = fastq.prepare_fastqs(
        str(r1),
        str(r2),
        None,
        str(fastp_dir),
        8,
    )

    assert processed_reads == fastq.FastqPaths(
        r1=f"{fastp_dir}/R1.fq.gz",
        r2=f"{fastp_dir}/R2.fq.gz",
        single=f"{fastp_dir}/RS.fq.gz",
    )
    assert calls[0]["unpaired"] == f"{fastp_dir}/unpaired.fq"
    assert calls[0]["cpus"] == 4

    with gzip.open(processed_reads.single, "rb") as infile:
        assert infile.read() == unpaired_content


def test_prepare_unpaired_reads_ignores_empty_unpaired_file(tmp_path):
    fastp_dir = tmp_path / "fastp"
    unpaired_file = fastp_dir / "unpaired.fq"

    fastp_dir.mkdir()
    unpaired_file.write_bytes(b"")

    assert fastq.prepare_unpaired_reads(str(unpaired_file), str(fastp_dir)) is None
    assert not (fastp_dir / "RS.fq").exists()
    assert not (fastp_dir / "RS.fq.gz").exists()
