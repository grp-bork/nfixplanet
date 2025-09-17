import logging
import sys

from . import pipeline


def main():
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        level=logging.DEBUG,
    )
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    logger.info("start pipeline")

    # TODO: parse args
    input = sys.argv[1]
    output = sys.argv[2]
    pipeline.run(input, output)


if __name__ == "__main__":
    main()
