import subprocess
import logging
import sys
from pathlib import Path

from nfixplanet.constants import (
    NFIXPLANET_CACHE_DIR,
    REFERENCE_FASTA_URL,
    REFERENCE_FASTA_NAME,
    REFERENCE_INDEX_NAME,
)

logger = logging.getLogger(__name__)


# ---------------------------
# Generic helpers
# ---------------------------


def check_files_exist(paths: list[str]):
    for p in paths:
        if not Path(p).exists():
            logger.error(f"Missing file: {p}")
            sys.exit(1)


def run_command(cmd: list[str], error_msg: str):
    try:
        logger.info("Running: " + " ".join(cmd))
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        logger.error(f"{error_msg}: {e}")
        sys.exit(1)


# ---------------------------
# FASTP
# ---------------------------


def fastp_paired(
    r1_in: str,
    r2_in: str,
    r1_out: str,
    r2_out: str,
    unpaired: str,
    html: str,
    json: str,
    cpus: int,
):
    run_command(
        [
            "fastp",
            "--thread",
            str(cpus),
            "--in1",
            r1_in,
            "--in2",
            r2_in,
            "--out1",
            r1_out,
            "--out2",
            r2_out,
            "--unpaired1",
            unpaired,
            "--unpaired2",
            unpaired,
            "--html",
            html,
            "--json",
            json,
            "--length_required",
            "45",
            "--cut_front",
            "--cut_tail",
            "--cut_window_size",
            "4",
            "--cut_mean_quality",
            "20",
        ],
        f"fastp paired failed on {r1_in} {r2_in}",
    )


def fastp_single(
    s_in: str,
    s_out: str,
    html: str,
    json: str,
    cpus: int,
):
    run_command(
        [
            "fastp",
            "--thread",
            str(cpus),
            "--in1",
            s_in,
            "--out1",
            s_out,
            "--html",
            html,
            "--json",
            json,
            "--length_required",
            "45",
            "--cut_front",
            "--cut_tail",
            "--cut_window_size",
            "4",
            "--cut_mean_quality",
            "20",
        ],
        f"fastp single failed on {s_in}",
    )


# ---------------------------
# HOSTILE
# ---------------------------


def hostile_index_fetch():
    run_command(
        ["hostile", "index", "fetch"],
        "hostile index fetch failed",
    )


def hostile_clean_paired(r1: str, r2: str, out_dir: str, cpus: int):
    run_command(
        [
            "hostile",
            "clean",
            "--fastq1",
            r1,
            "--fastq2",
            r2,
            "--threads",
            str(cpus),
            "--output",
            out_dir,
        ],
        f"hostile clean paired failed on {r1} {r2}",
    )

    r1_out = f"{out_dir}/R1.clean_1.fastq.gz"
    r2_out = f"{out_dir}/R2.clean_2.fastq.gz"
    return r1_out, r2_out


def hostile_clean_single(s, out_dir, cpus: int):
    run_command(
        [
            "hostile",
            "clean",
            "--fastq1",
            s,
            "--threads",
            str(cpus),
            "--output",
            out_dir,
        ],
        f"hostile clean single failed on {s}",
    )

    return f"{out_dir}/RS.clean.fastq.gz"


# ---------------------------
# MINIMAP2
# ---------------------------


def fetch_reference_file():
    # make sure cache dir exists
    NFIXPLANET_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    gz_path = f"{REFERENCE_FASTA_NAME}.gz"
    run_command(
        ["wget", "-O", gz_path, REFERENCE_FASTA_URL],
        "Failed to download reference file",
    )

    run_command(
        ["gunzip", gz_path],
        "Failed to unzip reference file",
    )


def minimap2_index(cpus: int):
    run_command(
        [
            "minimap2",
            "-d",
            str(REFERENCE_INDEX_NAME),
            "-t",
            str(cpus),
            str(REFERENCE_FASTA_NAME),
        ],
        "minimap2 index failed",
    )


def build_minimap_index(cpus: int):
    fetch_reference_file()
    # index
    minimap2_index(cpus)


# ---------------------------
# COVERM
# ---------------------------


def coverm_contig(
    r1: str | None,
    r2: str | None,
    single: str | None,
    reference_index: str,
    output_file: str,
    cpus: int,
):
    base_cmd = [
        "coverm",
        "contig",
        "-r",
        reference_index,
        "--threads",
        str(cpus),
        "--minimap2-reference-is-index",
        "--mapper",
        "minimap2-sr",
        "--methods",
        "covered_bases",
        "mean",
        "--output-file",
        output_file,
    ]

    if r1 and r2 and single:
        cmd = base_cmd + ["-1", r1, "-2", r2, "--single", single]
    elif r1 and r2:
        cmd = base_cmd + ["-1", r1, "-2", r2]
    elif single:
        cmd = base_cmd + ["--single", single]
    else:
        logger.error("No valid FASTQ inputs for coverm")
        sys.exit(1)

    run_command(cmd, "coverm failed")
