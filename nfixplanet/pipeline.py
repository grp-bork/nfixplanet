import sys
import logging

from . import external_tools
from . import processing

logger = logging.getLogger(__name__)


def run(input: str, output: str):
    external_tools.check_external_tools()

    # TODO: these need to go to temporary outputs
    #prodigal_output = "out/prodigal_out.fna"
    prodigal_output = "/scratch/robbani/repos/nfixplanet/reference_data/input/hmmscan/01_all_genes.fna"
    hmm_output = "out/hmm_out_all.tbl"

    # external_tools.prodigal(input, prodigal_output)
    #external_tools.hmmscan(prodigal_output, hmm_output)

    processing.filter(hmm_output)
