import os
import logging

from . import external_tools
from . import processing

logger = logging.getLogger(__name__)


def run(input: str, output_dir: str, genomic_context_range: int, cpus: int):
    external_tools.check_external_tools()

    prodigal_output = f"{output_dir}/prodigal_output.fna"
    hmm_output = f"{output_dir}/hmm_output.tbl"

    os.makedirs(output_dir, exist_ok=True)

    external_tools.prodigal(input, prodigal_output)
    external_tools.hmmscan(prodigal_output, hmm_output, cpus)

    processing.filter_and_write_files(hmm_output, output_dir, genomic_context_range)
    logger.info("Pipeline completed")
