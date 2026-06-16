import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path

from nfixplanet.profile import tools


@dataclass
class FastqPaths:
    r1: str | None
    r2: str | None
    single: str | None


def prepare_fastqs(
    r1: str | None,
    r2: str | None,
    single: str | None,
    fastp_dir: str,
    cpus: int,
) -> FastqPaths:
    processed_r1 = None
    processed_r2 = None
    processed_single = None
    unpaired_file = None

    if r1 and r2:
        processed_r1, processed_r2, unpaired_file = run_fastp_paired_reads(
            r1,
            r2,
            fastp_dir,
            cpus,
        )

    if single:
        singles_file = run_fastp_single_reads(single, fastp_dir, cpus)
        processed_single = merge_single_reads(singles_file, unpaired_file, fastp_dir)
    # TODO: handle unpaired_reads if no sinlge
    # elif unpaired_file:
    #     processed_single = unpaired_file

    return FastqPaths(processed_r1, processed_r2, processed_single)


def run_fastp_paired_reads(
    r1: str,
    r2: str,
    fastp_dir: str,
    cpus: int,
) -> tuple[str, str, str]:
    tools.check_files_exist([r1, r2])

    processed_r1 = f"{fastp_dir}/R1.fq.gz"
    processed_r2 = f"{fastp_dir}/R2.fq.gz"
    unpaired_file = f"{fastp_dir}/unpaired.fq"

    tools.fastp_paired(
        r1_in=r1,
        r2_in=r2,
        r1_out=processed_r1,
        r2_out=processed_r2,
        unpaired=unpaired_file,
        html=f"{fastp_dir}/paired.html",
        json=f"{fastp_dir}/paired.json",
        cpus=min(cpus, 4),
    )

    return processed_r1, processed_r2, unpaired_file


def run_fastp_single_reads(single: str, fastp_dir: str, cpus: int) -> str:
    tools.check_files_exist([single])

    singles_file = f"{fastp_dir}/singles.fq"
    tools.fastp_single(
        s_in=single,
        s_out=singles_file,
        html=f"{fastp_dir}/single.html",
        json=f"{fastp_dir}/single.json",
        cpus=cpus,
    )

    return singles_file


def merge_single_reads(
    singles_file: str,
    unpaired_file: str | None,
    fastp_dir: str,
) -> str | None:
    merged_reads = f"{fastp_dir}/RS.fq"
    merged_reads_gz = f"{fastp_dir}/RS.fq.gz"

    with open(merged_reads, "wb") as outfile:
        copy_file_if_present(singles_file, outfile)

        if unpaired_file:
            copy_file_if_present(unpaired_file, outfile)

    if not has_content(merged_reads):
        return None

    gzip_file(merged_reads, merged_reads_gz)
    return merged_reads_gz


def copy_file_if_present(path: str, outfile):
    if Path(path).exists():
        with open(path, "rb") as infile:
            shutil.copyfileobj(infile, outfile)


def has_content(path: str) -> bool:
    file_path = Path(path)
    return file_path.exists() and file_path.stat().st_size > 0


def gzip_file(input_path: str, output_path: str):
    with open(input_path, "rb") as f_in, gzip.open(output_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
