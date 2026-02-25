import os
import logging
import shutil

logger = logging.getLogger(__name__)


def configure_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    format = (
        "%(asctime)s : [%(levelname)7s] : %(name)s:%(lineno)s %(funcName)20s() : %(message)s"
        if verbose
        else "[%(asctime)s] %(levelname)s: %(message)s"
    )

    logging.basicConfig(
        format=format,
        datefmt="%H:%M:%S",
        level=level,
        force=True,
    )


def check_files_exist(files: list[str]) -> None:
    for file in files:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"{file} does not exist")


def check_if_tool_exists(tool_name: str) -> bool:
    """Check if tool is available."""
    return shutil.which(tool_name) is not None


def check_external_tools(required_tools: list[str]):
    """Checks if tool dependencies are available."""
    for tool in required_tools:
        if not check_if_tool_exists(tool):
            logger.error(f"{tool} not found.")
            exit(1)

    logger.info("All required software found")
