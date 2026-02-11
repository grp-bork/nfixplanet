import os
import pandas as pd
import pytest

from nfixplanet import processing
from nfixplanet.constants import GENE_FAMILIES


@pytest.mark.parametrize("gene_family", GENE_FAMILIES, ids=lambda f: f.name)
def test_filter_top_hits_by_genes(gene_family):
    # load hmm_output once from a known test file
    hmm_path = os.path.join("tests", "references", "input", "filter", "hmm_out_all.tbl")
    hmm_output = processing.load_hmm_output(hmm_path)
    best_hits = processing.get_best_hmm_hits(hmm_output)

    # run filtering
    genes_to_hits = processing.filter_top_hits_by_genes(best_hits, gene_family, genomic_context_range=10)

    assert genes_to_hits is not None, "All test data should return hits"

    for gene, df in genes_to_hits.items():
        if df.empty:
            continue

        if gene_family.name in ("nifE", "nifN") and gene in ("nifD", "nifH", "nifK"):
            continue

        # drop extra cols for fair comparison
        df = df.drop(["contig", "gene"], axis=1)

        # path to reference
        ref_path = os.path.join("tests", "references", "output", "filter", f"{gene}.tsv")
        assert os.path.exists(ref_path), f"Missing reference file: {ref_path}"

        # read reference, no header
        ref_df = pd.read_csv(ref_path, sep="\t", header=None)

        # reset index for clean comparison
        df_values = df.reset_index(drop=True).values
        ref_values = ref_df.reset_index(drop=True).values

        assert (
            df_values == ref_values
        ).all(), f"Mismatch for {gene_family.name}/{gene}"
        break
