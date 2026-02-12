# install with pip install -e .

# nfixplanet --input_genomes tests/references/input/prodigal/GCA_001049335.1.genomes_clean.fa --output_directory out/test3 --verbose
# nfixplanet --input_orfs tests/references/input/hmmscan/GCA_001049335.1.genomes_clean.fna --output_directory out/test4 --verbose
nfixplanet --input_hmms tests/references/input/filter/hmm_out_all.tbl --output_directory out/test5 --verbose
