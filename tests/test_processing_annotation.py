import os
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from nfixplanet.annotate import processing
from nfixplanet.constants import GENE_FAMILIES


@pytest.mark.parametrize("gene_family", GENE_FAMILIES, ids=lambda f: f.name)
def test_filter_top_hits(gene_family):
    # load hmm_output once from a known test file
    hmm_path = os.path.join("tests", "references", "input", "filter", "hmm_output.tbl.gz")
    hmm_output = processing.load_hmm_output(hmm_path)
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
