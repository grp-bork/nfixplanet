import argparse
import logging
import sys
from pathlib import Path

from . import pipeline


def parse_args(argv: list[str]):
    parser = argparse.ArgumentParser(
        description="A pipeline for the detection of nitrogen fixers.\n"
    )

    parser.add_argument(
        "--input_fasta",
        help="Path to input fasta file (can be gzipped)",
    )

    parser.add_argument(
        "--output_directory",
        help="Path to output directory",
    )

    # resolve default hmm profile inside package
    default_hmm = (
        Path(__file__).parent.parent
        / "reference_data"
        / "hmm_profiles"
        / "nfixplanet_models.hmm"
    )
    parser.add_argument(
        "--hmm_profile_path",
        default=default_hmm,
        help=f"Path to HMM profile file (default: {default_hmm})",
    )

    parser.add_argument(
        "--cpus",
        type=int,
        default=2,
        help="Number of CPUs to use for HMMscan (default: 2, max recommended: 4)",
    )

    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )

    parser.add_argument(
        "--version",
        action="version",
        help="Print version number and exit.",
        version="0.1.0", # TODO: dynamically update version
    )

    return parser.parse_args(argv)


def main():
    args = parse_args(sys.argv[1:])
    if args.verbose:
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
    logger = logging.getLogger(__name__)
    logger.debug(args)
    logger.info("Start pipeline")
    pipeline.run(
        args.input_fasta, args.output_directory, args.hmm_profile_path, args.cpus
    )


if __name__ == "__main__":
    main()
