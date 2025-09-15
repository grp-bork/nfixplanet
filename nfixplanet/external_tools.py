import shutil
import subprocess
import logging
import sys

logger = logging.getLogger(__name__)

from .config import HMM_PROFILE_PATH, CPUS


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


def hmmscan(input_file: str, out_file: str):
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
                str(CPUS),
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


def check_if_tool_exists(tool_name: str) -> bool:
    """Check if tool is available."""
    return shutil.which(tool_name) is not None
