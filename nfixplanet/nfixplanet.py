import logging
import sys

from . import pipeline

def main():

    logger = logging.getLogger("nfixplanet.py")
    # TODO: parse args
    input = sys.argv[1]
    output = sys.argv[2]
    pipeline.run_pipeline(input, output)

if __name__ == "__main__":
    main()
