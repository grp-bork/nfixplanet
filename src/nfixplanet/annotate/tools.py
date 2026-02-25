import subprocess
import logging
import sys

from nfixplanet.constants import HMM_PROFILE_PATH

logger = logging.getLogger(__name__)


def prodigal(input_file: str, out_file: str):
    """Run prodigal

    Arguments:
        input_file (str): full path of input fasta
        out_file (str): fullpath of output file
    """
    logger.info("Running prodigal...")
    try:
        if input_file.endswith(".gz"):
            uncompressed = subprocess.Popen(
                ("zcat", input_file), stdout=subprocess.PIPE
            )
            subprocess.check_output(
                [
                    "prodigal",
                    "-i",
                    "/dev/stdin",
                    "-d",
                    out_file,
                    "-o",
                    "/dev/null",
                    "-p",
                    "meta",
                    "-q",
                ],
                stdin=uncompressed.stdout,
                universal_newlines=True,
            )
        else:
            subprocess.check_output(
                [
                    "prodigal",
                    "-i",
                    input_file,
                    "-d",
                    out_file,
                    "-o",
                    "/dev/null",
                    "-p",
                    "meta",
                    "-q",
                ],
                universal_newlines=True,
            )
        logger.info("Prodigal completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run Prodigal on {input_file}: {e}")
        sys.exit(1)


def hmmscan(input_file: str, out_file: str, cpus: int):
    """Run hmmscan

    Arguments:
        input_file (str): full path of input fasta (fna)
        out_file (str): fullpath of output file
    """
    logger.info("Running hmmscan...")
    try:
        subprocess.check_output(
            [
                "hmmscan",
                "--cpu",
                str(cpus),
                "--cut_ga",
                "--tblout",
                out_file,
                HMM_PROFILE_PATH,
                input_file,
            ],
            universal_newlines=True,
        )
        logger.info("Hmmscan completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run hmmscan on {input_file}: {e}")
        sys.exit(1)


def fastp_paired(r1_in: str, r2_in: str, r1_out: str, r2_out: str, cpus: int):
    logger.info("Running fastp (paired)...")
    try:
        subprocess.check_output(
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
                "--length_required 45",
                "--cut_front",
                "--cut_tail",
                "--cut_window_size",
                "4",
                "--cut_mean_quality",
                "20",
            ],
            universal_newlines=True,
        )
        logger.info("Hmmscan completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run hmmscan on {r1_in} & {r2_in}: {e}")
        sys.exit(1)


def fastp_single(s_in: str, s_out: str, cpus: int):
    logger.info("Running fastp (single)...")
    try:
        subprocess.check_output(
            [
                "fastp",
                "--thread",
                str(cpus),
                "--in1",
                s_in,
                "--out1",
                s_out,
                "--length_required 45",
                "--cut_front",
                "--cut_tail",
                "--cut_window_size",
                "4",
                "--cut_mean_quality",
                "20",
            ],
            universal_newlines=True,
        )
        logger.info("Hmmscan completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run hmmscan on {s_in} & {s_out}: {e}")
        sys.exit(1)
