import os
import filecmp
import pytest
from pathlib import Path

from nfixplanet import external_tools


@pytest.mark.parametrize(
    "input_file,ref_file",
    [
        (
            os.path.join(
                "tests",
                "references",
                "input",
                "prodigal",
                "GCA_001049335.1.genomes_clean.fa",
            ),
            os.path.join(
                "tests",
                "references",
                "output",
                "prodigal",
                "GCA_001049335.1.genomes.genes.fna",
            ),
        ),
    ],
)
def test_prodigal(tmp_path, input_file, ref_file):
    output_file = tmp_path / "prodigal_out.fna"

    external_tools.prodigal(str(input_file), str(output_file))

    assert os.path.exists(output_file), "Prodigal did not create output file"
    assert os.path.exists(ref_file), f"Missing reference file: {ref_file}"

    # compare bodies
    assert filecmp.cmp(output_file, ref_file, shallow=False), "Prodigal output mismatch"


@pytest.mark.parametrize(
    "input_file,ref_file",
    [
        (
            os.path.join(
                "tests",
                "references",
                "input",
                "hmmscan",
                "GCA_001049335.1.genomes_clean.fna",
            ),
            os.path.join(
                "tests",
                "references",
                "output",
                "hmmscan",
                "GCA_001049335.1.genomes.genes.tbl",
            ),
        ),
    ],
)
def test_hmmscan(tmp_path, input_file, ref_file):
    output_file = tmp_path / "hmm_out_all.tbl"
    hmm_profile_path = (
        Path(__file__).parent.parent
        / "reference_data"
        / "hmm_profiles"
        / "nfixplanet_models.hmm"
    )

    external_tools.hmmscan(str(input_file), str(output_file), str(hmm_profile_path), 2)

    assert os.path.exists(output_file), "hmmscan did not create output file"
    assert os.path.exists(ref_file), f"Missing reference file: {ref_file}"

    def strip_comments(path):
        with open(path) as f:
            return [line for line in f.read().splitlines() if not line.startswith("#")]

    out_lines = strip_comments(output_file)
    ref_lines = strip_comments(ref_file)

    assert out_lines == ref_lines, "hmmscan output mismatch"
