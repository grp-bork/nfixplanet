import os
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from nfixplanet.annotate import processing
from nfixplanet.constants import GENE_FAMILIES


def test_load_hmm_output_custom_format_swaps_target_and_query():
    hmm_path = os.path.join(
        "tests", "references", "input", "filter_custom", "nfix_issue2919_clean.tsv"
    )

    hmm_output = processing.load_hmm_output(hmm_path, custom_hmm_format=True)

    assert "inc" in hmm_output.columns
    assert hmm_output.loc[0, "query_name"] == (
        "SAMEA112550948.psa_megahit.prodigal.fna;k141_2452185_1"
    )
    assert hmm_output.loc[0, "target_name"] == "vnfH"
    assert hmm_output.loc[0, "contig"] == (
        "SAMEA112550948.psa_megahit.prodigal.fna;k141_2452185"
    )
    assert hmm_output.loc[0, "gene"] == 1


def test_filter_custom_hmm_output_runs():
    hmm_path = os.path.join(
        "tests", "references", "input", "filter_custom", "nfix_issue2919_clean.tsv"
    )
    hmm_output = processing.load_hmm_output(hmm_path, custom_hmm_format=True)
    best_hits = processing.get_best_hmm_hits(hmm_output)

    non_empty_hits = []
    for gene_family in GENE_FAMILIES:
        genes_to_hits = processing.get_filtered_top_hits(
            best_hits, gene_family, genomic_context_range=10
        )
        if genes_to_hits is None:
            continue

        non_empty_hits.extend(
            gene_df for gene_df in genes_to_hits.values() if not gene_df.empty
        )

    assert non_empty_hits, "Custom HMM input should produce at least one filtered hit"


@pytest.mark.parametrize("gene_family", GENE_FAMILIES, ids=lambda f: f.name)
def test_filter_top_hits(gene_family):
    # load hmm_output once from a known test file
    hmm_path = os.path.join(
        "tests", "references", "input", "filter", "hmm_output.tbl.gz"
    )
    hmm_output = processing.load_hmm_output(hmm_path, custom_hmm_format=False)
    best_hits = processing.get_best_hmm_hits(hmm_output)

    # run filtering
    genes_to_hits = processing.get_filtered_top_hits(
        best_hits, gene_family, genomic_context_range=10
    )

    assert genes_to_hits is not None, "All test data should return hits"

    for gene, gen_df in genes_to_hits.items():
        if gen_df.empty:
            continue

        # sort for fair comparison
        gen_df = gen_df.sort_values(by=gen_df.columns[0]).reset_index(drop=True)

        # path to reference
        ref_path = os.path.join(
            "tests", "references", "output", "filter", f"{gene}.tsv"
        )
        assert os.path.exists(ref_path), f"Missing reference file: {ref_path}"

        # read reference, no header
        ref_df = pd.read_csv(
            ref_path,
            sep="\t",
        )
        ref_df = ref_df.sort_values(by=ref_df.columns[0]).reset_index(drop=True)

        # print(gene)
        # print(ref_df)
        # print(gen_df)
        assert_frame_equal(ref_df, gen_df, check_dtype=False)

        # try:
        #     assert_frame_equal(ref_df, gen_df, check_dtype=False)
        # except AssertionError:
        #     print(f"failed: {gene}")
        #     gen_df.to_csv(f"{gene}.tsv", sep="\t", index=False)
