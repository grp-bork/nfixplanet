import logging
from pathlib import Path

from nfixplanet.profile import fastq
from nfixplanet.profile import tools
from nfixplanet.profile import processing
from nfixplanet import utils
from nfixplanet.constants import (
    MAP_TOOLS,
    HOSTILE_CACHE_DIR,
    NFIXPLANET_CACHE_DIR,
    REFERENCE_INDEX_NAME,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)


def prepare_work_dirs(work_dir: str) -> tuple[str, str]:
    ensure_dir(work_dir)

    fastp_dir = f"{work_dir}/processed_fastqs"
    clean_dir = f"{work_dir}/cleaned_fastqs"

    ensure_dir(fastp_dir)
    ensure_dir(clean_dir)

    return fastp_dir, clean_dir


def run_profile(
    sample_id: str,
    r1: str | None,
    r2: str | None,
    single: str | None,
    output_dir: str,
    work_dir: str,
    cpus: int,
    input_coverage: str | None,
):
    """
    preprocess_fastqs -> clean_fastq -> coverm -> process results
    OR coverm -> process results
    """
    ensure_dir(output_dir)

    if input_coverage:
        tools.check_files_exist([input_coverage])
        logger.info(
            "Input CoverM coverage provided, skipping FASTQ processing and mapping"
        )
        process_coverage_results(sample_id, input_coverage, output_dir)
        logger.info(f"Finished processing coverage table for sample {sample_id}")
        return

    utils.check_external_tools(MAP_TOOLS)

    # Fastp
    fastp_dir, clean_dir = prepare_work_dirs(work_dir)
    processed_reads = fastq.prepare_fastqs(r1, r2, single, fastp_dir, cpus)

    # Hostile
    ensure_hostile_index()
    cleaned_reads = clean_fastqs(processed_reads, clean_dir, cpus)
    
    # CoverM
    ensure_reference_index(cpus)
    coverage_file = run_coverm(sample_id, output_dir, cleaned_reads, cpus)
    
    logger.info(f"Finished calculating coverage score for sample {sample_id}")

    # Processing
    process_coverage_results(sample_id, coverage_file, output_dir)

    logger.info("Pipeline complete")


def ensure_hostile_index():
    hostile_cache_path = Path(HOSTILE_CACHE_DIR).expanduser()
    hostile_cache_path.mkdir(parents=True, exist_ok=True)
    if not any(hostile_cache_path.iterdir()):
        tools.hostile_index_fetch()


def clean_fastqs(
    processed_reads: fastq.FastqPaths,
    clean_dir: str,
    cpus: int,
) -> fastq.FastqPaths:
    cleaned_r1 = None
    cleaned_r2 = None
    cleaned_single = None

    if processed_reads.r1 and processed_reads.r2:
        cleaned_r1, cleaned_r2 = tools.hostile_clean_paired(
            r1=processed_reads.r1,
            r2=processed_reads.r2,
            out_dir=clean_dir,
            cpus=cpus,
        )

    if processed_reads.single:
        cleaned_single = tools.hostile_clean_single(
            s=processed_reads.single,
            out_dir=clean_dir,
            cpus=cpus,
        )

    return fastq.FastqPaths(cleaned_r1, cleaned_r2, cleaned_single)


def ensure_reference_index(cpus: int):
    NFIXPLANET_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    if not REFERENCE_INDEX_NAME.exists():
        tools.build_minimap_index(min(cpus, 3))


def run_coverm(
    sample_id: str,
    output_dir: str,
    cleaned_reads: fastq.FastqPaths,
    cpus: int,
) -> str:
    coverage_file = f"{output_dir}/{sample_id}_coverage.tsv"

    tools.coverm_contig(
        r1=cleaned_reads.r1,
        r2=cleaned_reads.r2,
        single=cleaned_reads.single,
        reference_index=str(REFERENCE_INDEX_NAME),
        output_file=coverage_file,
        cpus=cpus,
    )

    return coverage_file


def process_coverage_results(sample_id: str, coverage_file: str, output_dir: str):
    logger.info("Processing CoverM results")
    gene_file = f"{output_dir}/{sample_id}_gene_table.tsv"
    otu_file = f"{output_dir}/{sample_id}_OTU_table.tsv"
    processing.annotate_mapping_results(
        coverage_file,
        gene_file,
        otu_file,
        sample_name=sample_id,
    )
    processing.profile_results(otu_file, output_dir)
