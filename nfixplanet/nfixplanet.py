import logging
import sys

from . import external_tools

logger = logging.getLogger("nfixplanet.py")

def start_checks():
    """Checks if tool dependencies are available."""
    if not external_tools.check_if_tool_exists("prodigal"):
        logger.error("Prodigal not found.")
        sys.exit(1)
    # if not external_tools.check_if_tool_exists("grep"):
    #     logger.error("grep not found.")
    #     sys.exit(1)
    # if not external_tools.check_if_tool_exists("zcat"):
    #     logger.error("zcat not found.")
    #     sys.exit(1)

def main():
    # TODO: parse args
    start_checks()

if __name__ == "__main__":
    main()
