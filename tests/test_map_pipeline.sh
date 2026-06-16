#!/bin/bash

# This test just checks that the CoverM coverage files can be produced
# All other tests are in the `test_profile_pipeline.py` file


mkdir -p tmp/work
mkdir -p tmp/output

nfixplanet map --read_1 tests/references/input/map_pipeline/SRR5371433.pair.1.fq.gz \
               --read_2 tests/references/input/map_pipeline/SRR5371433.pair.2.fq.gz \
               --single tests/references/input/map_pipeline/SRR5371433.single.fq.gz \
               --sample_id SAMN06117828 \
               --work_directory tmp/work \
               --output_directory tmp/output
diff -s tmp/output/SAMN06117828_coverage.tsv tests/references/output/map_pipeline/SAMN06117828_coverage.tsv

# rm -r tmp
