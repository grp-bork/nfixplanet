import os
import logging

from nfixplanet.annotate import tools
from nfixplanet.annotate import processing
from nfixplanet import utils
from nfixplanet.constants import ANNOTATE_TOOLS

logger = logging.getLogger(__name__)


def run_annotate(
    input_genome: str,
    input_orf: str,
    input_hmm: str,
    output_dir: str,
    genomic_context_range: int,
    cpus: int,
):
    """Run nitrogen fixer annotation pipeline"""
    utils.check_external_tools(ANNOTATE_TOOLS)
    os.makedirs(output_dir, exist_ok=True)

    if input_genome:
        genome_path = input_genome
        orf_path = f"{output_dir}/prodigal_output.fna"
        hmm_path = f"{output_dir}/hmm_output.tbl"

        tools.prodigal(genome_path, orf_path)
        tools.hmmscan(orf_path, hmm_path, cpus)
    elif input_orf:
        orf_path = input_orf
        hmm_path = f"{output_dir}/hmm_output.tbl"
        logger.info("Input ORFs provided, skipping Prodigal")

        tools.hmmscan(orf_path, hmm_path, cpus)
    elif input_hmm:
        hmm_path = input_hmm
        logger.info("Input HMMs provided, skipping Prodigal and HMMscan")

    processing.filter_and_write_files(hmm_path, output_dir, genomic_context_range)
    logger.info("Pipeline completed")
