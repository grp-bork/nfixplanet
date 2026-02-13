import os
import logging


def configure_logging(verbose: bool):
    if verbose:
        logging.basicConfig(
            format="%(asctime)s : [%(levelname)7s] : %(name)s:%(lineno)s %(funcName)20s() : %(message)s",
            datefmt="%H:%M:%S",
            level=logging.DEBUG,
        )
    else:
        logging.basicConfig(
            format="[%(asctime)s] %(levelname)s: %(message)s",
            datefmt="%H:%M:%S",
            level=logging.INFO,
        )


def check_files_exist(files: list[str]) -> None:
    for file in files:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"{file} does not exist")
