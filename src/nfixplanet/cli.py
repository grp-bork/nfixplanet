import argparse
import logging
import sys
import tomllib
from pathlib import Path

from . import pipeline


def get_version() -> str:
    pyproject_path = Path(__file__).parent.parent.parent / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        pyproject_data = tomllib.load(f)
    return pyproject_data["project"]["version"]


def get_parser():
    parser = argparse.ArgumentParser(
        description="A pipeline for the detection of nitrogen fixers.\n"
    )

    parser.add_argument(
        "--input_genomes",
        help="Path to input genom fasta file (can be gzipped)",
        default="",
    )

    parser.add_argument(
        "--input_orfs",
        help="Path to input ORF fasta file (Prodigal output)",
        default="",
    )

    parser.add_argument(
        "--input_hmms",
        help="Path to HMM tables (HMMER output)",
        default="",
    )

    parser.add_argument(
        "--output_directory",
        help="Path to output directory",
        required=True,
    )

    parser.add_argument(
        "--genomic_context_range",
        type=int,
        default=10,
        help="Maximum number of genes upstream or downstream to consider for the operon context (default: 10)",
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
        version=get_version(),
    )

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args(sys.argv[1:])
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

    if not (args.input_genomes or args.input_orfs or args.input_hmms):
        parser.error(
            "At least one of --input_genomes, --input_orfs, or --input_hmms must be provided"
        )

    logger.info("Start pipeline")
    pipeline.run(
        args.input_genomes,
        args.input_orfs,
        args.input_hmms,
        args.output_directory,
        args.genomic_context_range,
        args.cpus,
    )


if __name__ == "__main__":
    main()
