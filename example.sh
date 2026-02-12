# install with pip install -e .

# nfixplanet --input_genomes tests/references/input/prodigal/GCA_001049335.1.genomes_clean.fa --output_directory out/full --verbose
# nfixplanet --input_orfs tests/references/input/hmmscan/GCA_001049335.1.genomes_clean.fna --output_directory out/cut_ga --verbose
# nfixplanet --input_orfs tests/references/input/hmmscan/GCA_001049335.1.genomes_clean.fna --output_directory out/no_ga --verbose
# nfixplanet --input_hmms tests/references/input/filter/hmm_out_all.tbl --output_directory out/filter --verbose
