import logging
from pathlib import Path

from . import tools

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)


def run_map(
    sample_id: str,
    r1: str | None,
    r2: str | None,
    single: str | None,
    reference_index: str,
    work_dir: str | None = None,
    cpus: int = 4,
):
    """
    preprocess_fastqs -> clean_fastq -> coverm
    """
    #TODO: if not reference/index build them

    tmp_dir = work_dir if work_dir else f"tmp_{sample_id}"
    ensure_dir(tmp_dir)

    fastp_dir = f"{tmp_dir}/processed_fastqs"
    clean_dir = f"{tmp_dir}/cleaned_fastqs"

    ensure_dir(fastp_dir)
    ensure_dir(clean_dir)

    processed_r1 = None
    processed_r2 = None
    processed_s = None

    # --------------------
    # FASTP
    # --------------------
    if r1 and r2:
        tools.check_files_exist([r1, r2])
        processed_r1 = f"{fastp_dir}/R1.fq.gz"
        processed_r2 = f"{fastp_dir}/R2.fq.gz"

        tools.fastp_paired(
            r1_in=r1,
            r2_in=r2,
            r1_out=processed_r1,
            r2_out=processed_r2,
            cpus=cpus,
        )

    if single:
        tools.check_files_exist([single])
        processed_s = f"{fastp_dir}/RS.fq.gz"

        tools.fastp_single(
            s_in=single,
            s_out=processed_s,
            cpus=cpus,
        )

    # --------------------
    # HOSTILE CLEAN
    # --------------------
    cleaned_r1 = None
    cleaned_r2 = None
    cleaned_s = None

    if processed_r1 and processed_r2:
        cleaned_r1, cleaned_r2 = tools.hostile_clean_paired(
            processed_r1,
            processed_r2,
            clean_dir,
            cpus=cpus,
        )

    if processed_s:
        cleaned_s = tools.hostile_clean_single(
            processed_s,
            clean_dir,
            cpus=cpus,
        )

    # --------------------
    # COVERM
    # --------------------
    output_file = f"{sample_id}_sample_coverage.tsv"

    tools.coverm_contig(
        r1=cleaned_r1,
        r2=cleaned_r2,
        single=cleaned_s,
        reference_index=reference_index,
        output_file=output_file,
        cpus=cpus,
    )

    logger.info(f"Finished sample {sample_id}")
    return output_file


def build_minimap_index(reference_fasta: str, cpus: int = 4) -> str:
    # TODO: provide option for reference file
    tools.check_files_exist([reference_fasta])
    index_path = f"{Path(reference_fasta).stem}.mmi"
    tools.minimap2_index(reference_fasta, index_path, cpus)
    return index_path


def build_hostile_index():
    # TODO: provide option for index?
    tools.hostile_index_fetch()
