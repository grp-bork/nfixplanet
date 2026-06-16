import gzip
import shutil
import logging
from pathlib import Path

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
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if input_coverage:
        tools.check_files_exist([input_coverage])
        logger.info(
            "Input CoverM coverage provided, skipping FASTQ processing and mapping"
        )
        process_coverage_results(sample_id, input_coverage, output_dir)
        logger.info(f"Finished processing coverage table for sample {sample_id}")
        return

    utils.check_external_tools(MAP_TOOLS)
    ensure_dir(work_dir)

    fastp_dir = f"{work_dir}/processed_fastqs"
    clean_dir = f"{work_dir}/cleaned_fastqs"

    ensure_dir(fastp_dir)
    ensure_dir(clean_dir)

    processed_r1 = None
    processed_r2 = None
    processed_s = None

    # --------------------
    # FASTP
    # --------------------
    unpaired_file = None
    singles_file = None

    if r1 and r2:
        tools.check_files_exist([r1, r2])

        processed_r1 = f"{fastp_dir}/R1.fq.gz"
        processed_r2 = f"{fastp_dir}/R2.fq.gz"
        unpaired_file = f"{fastp_dir}/unpaired.fq"

        paired_html = f"{fastp_dir}/paired.html"
        paired_json = f"{fastp_dir}/paired.json"

        tools.fastp_paired(
            r1_in=r1,
            r2_in=r2,
            r1_out=processed_r1,
            r2_out=processed_r2,
            unpaired=unpaired_file,
            html=paired_html,
            json=paired_json,
            cpus=min(cpus, 4),
        )

    if single:
        tools.check_files_exist([single])

        singles_file = f"{fastp_dir}/singles.fq"

        single_html = f"{fastp_dir}/single.html"
        single_json = f"{fastp_dir}/single.json"

        tools.fastp_single(
            s_in=single,
            s_out=singles_file,
            html=single_html,
            json=single_json,
            cpus=cpus,
        )
        # --------------------
        # MERGE UNPAIRED + SINGLES
        # --------------------
        merged_rs = f"{fastp_dir}/RS.fq"
        merged_rs_gz = f"{fastp_dir}/RS.fq.gz"

        with open(merged_rs, "wb") as outfile:
            if singles_file and Path(singles_file).exists():
                with open(singles_file, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

            if unpaired_file and Path(unpaired_file).exists():
                with open(unpaired_file, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

        # gzip only if file has content
        if Path(merged_rs).exists() and Path(merged_rs).stat().st_size > 0:
            with open(merged_rs, "rb") as f_in, gzip.open(merged_rs_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

            processed_s = merged_rs_gz
        else:
            processed_s = None

    # --------------------
    # HOSTILE INDEX
    # --------------------
    hostile_cache_path = Path(HOSTILE_CACHE_DIR).expanduser()
    hostile_cache_path.mkdir(parents=True, exist_ok=True)
    if not any(hostile_cache_path.iterdir()):
        tools.hostile_index_fetch()

    # --------------------
    # HOSTILE CLEAN
    # --------------------
    cleaned_r1 = None
    cleaned_r2 = None
    cleaned_s = None

    if processed_r1 and processed_r2:
        cleaned_r1, cleaned_r2 = tools.hostile_clean_paired(
            r1=processed_r1,
            r2=processed_r2,
            out_dir=clean_dir,
            cpus=cpus,
        )

    if processed_s:
        cleaned_s = tools.hostile_clean_single(
            s=processed_s,
            out_dir=clean_dir,
            cpus=cpus,
        )

    # --------------------
    # MINIMAP INDEX
    # --------------------
    NFIXPLANET_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    if not REFERENCE_INDEX_NAME.exists():
        tools.build_minimap_index(min(cpus, 3))

    # --------------------
    # COVERM
    # --------------------
    coverage_file = f"{output_dir}/{sample_id}_coverage.tsv"

    tools.coverm_contig(
        r1=cleaned_r1,
        r2=cleaned_r2,
        single=cleaned_s,
        reference_index=str(REFERENCE_INDEX_NAME),
        output_file=coverage_file,
        cpus=cpus,
    )

    process_coverage_results(sample_id, coverage_file, output_dir)

    logger.info(f"Finished calculating coverage score for sample {sample_id}")


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
