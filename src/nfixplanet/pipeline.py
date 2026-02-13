import os
import logging

from . import external_tools
from . import processing
from . import utils

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
    external_tools.check_external_tools()
    os.makedirs(output_dir, exist_ok=True)

    if input_genome:
        genome_path = input_genome
        orf_path = f"{output_dir}/prodigal_output.fna"
        hmm_path = f"{output_dir}/hmm_output.tbl"

        external_tools.prodigal(genome_path, orf_path)
        external_tools.hmmscan(orf_path, hmm_path, cpus)
    elif input_orf:
        orf_path = input_orf
        hmm_path = f"{output_dir}/hmm_output.tbl"
        logger.info("Input ORFs provided, skipping Prodigal")

        external_tools.hmmscan(orf_path, hmm_path, cpus)
    elif input_hmm:
        hmm_path = input_hmm
        logger.info("Input HMMs provided, skipping Prodigal and HMMscan")

    processing.filter_and_write_files(hmm_path, output_dir, genomic_context_range)
    logger.info("Pipeline completed")


def run_map(r1: str, r2: str, s: str, work_dir: str | None = None):
    tmp_dir = work_dir if work_dir else "tmp"
    os.mkdir(tmp_dir)

    if r1 and r2:
        utils.check_files_exist([r1, r2])
        external_tools.fastp_paired(r1, r2, f"{tmp_dir}/r1.fq", f"{tmp_dir}/r2.fq")

    if s:
        pass
    # run fastp
    # run hostile
    # hostile index
    # hostile clean
    # run minimap
    # minimap index
    # run coverm
    # run lucas R scripts
    pass
