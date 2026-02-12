import argparse
import logging
import sys
import tomllib
from pathlib import Path

from . import pipeline, utils


def get_version() -> str:
    pyproject_path = Path(__file__).parent.parent.parent / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        pyproject_data = tomllib.load(f)
    return pyproject_data["project"]["version"]


def get_parser():
    parser = argparse.ArgumentParser(
        description="A pipeline for the detection of nitrogen fixers."
    )

    parser.add_argument(
        "--version",
        action="version",
        help="Print version number and exit.",
        version=get_version(),
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # -------------------
    # annotate subcommand
    # -------------------
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Run annotation pipeline",
    )

    annotate_parser.add_argument(
        "--input_genomes",
        help="Path to input genom fasta file (can be gzipped)",
        default="",
    )
    annotate_parser.add_argument(
        "--input_orfs",
        help="Path to input ORF fasta file (Prodigal output)",
        default="",
    )

    annotate_parser.add_argument(
        "--input_hmms",
        help="Path to HMM tables (HMMER output)",
        default="",
    )

    annotate_parser.add_argument(
        "--output_directory",
        required=True,
        help="Path to output directory",
    )

    annotate_parser.add_argument(
        "--genomic_context_range",
        type=int,
        default=10,
        help="Maximum number of genes upstream or downstream to consider for the operon context (default: 10)",
    )

    annotate_parser.add_argument(
        "--cpus",
        type=int,
        default=2,
        help="Number of CPUs to use for HMMscan (default: 2, max recommended: 4)",
    )

    annotate_parser.add_argument(
        "--verbose",
        action="store_true",
    )

    annotate_parser.set_defaults(func=run_annotate_command)

    # ---------------
    # map subcommand
    # ---------------
    map_parser = subparsers.add_parser(
        "map",
        help="Run mapping pipeline",
    )

    map_parser.add_argument(
        "--output_directory",
        required=True,
    )

    map_parser.add_argument(
        "--verbose",
        action="store_true",
    )

    map_parser.set_defaults(func=run_map_command)

    return parser


def run_annotate_command(args):
    utils.configure_logging(args.verbose)
    logger = logging.getLogger(__name__)
    logger.debug(args)

    if not (args.input_genomes or args.input_orfs or args.input_hmms):
        raise SystemExit(
            "At least one of --input_genomes, --input_orfs, or --input_hmms must be provided"
        )

    logger.info("Start annotation pipeline")

    pipeline.run_annotate(
        args.input_genomes,
        args.input_orfs,
        args.input_hmms,
        args.output_directory,
        args.genomic_context_range,
        args.cpus,
    )


def run_map_command(args):
    utils.configure_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Start mapping pipeline")

    pipeline.run_map()


def main():
    parser = get_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
