#!/bin/bash

nfixplanet map --read_1 tests/references/input/map_pipeline/SRR5371433.pair.1.fq.gz \
               --read_2 tests/references/input/map_pipeline/SRR5371433.pair.2.fq.gz \
               --single tests/references/input/map_pipeline/SRR5371433.single.fq.gz \
               --sample_id SAMN06117828 \
               --work_directory tmp \
               --output_directory tmp
diff -s tmp/SAMN06117828_sample_coverage.tsv tests/references/output/map_pipeline/SAMN06117828_sample_coverage.tsv
