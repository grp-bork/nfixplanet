import sys
import logging

from . import external_tools

logger = logging.getLogger(__name__)

def start_checks():
    """Checks if tool dependencies are available."""
    if not external_tools.check_if_tool_exists("prodigal"):
        logger.error("Prodigal not found.")
        sys.exit(1)
    if not external_tools.check_if_tool_exists("hmmscan"):
        logger.error("Hmmscan not found.")
        sys.exit(1)
    # if not external_tools.check_if_tool_exists("grep"):
    #     logger.error("grep not found.")
    #     sys.exit(1)
    if not external_tools.check_if_tool_exists("zcat"):
        logger.error("zcat not found.")
        sys.exit(1)
    logger.info("All required software found")

def run_pipeline(input:str, output:str):
    start_checks()

    # TODO: these need to go to temporary outputs
    prodigal_output = "out/prodigal_out.fna"
    hmm_output = "out/hmm_out.tbl"
    
    logger.info("Running prodigal...")
    external_tools.prodigal(input, prodigal_output)
    logger.info("Prodigal completed")
    logger.info("Running hmmscan...")
    external_tools.hmmscan(prodigal_output, hmm_output)
    logger.info("Hmmscan completed")
