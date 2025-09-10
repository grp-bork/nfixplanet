#!/bin/bash

# run using ./01_annotation_pipeline.sh
# activate conda env before running
#conda activate nfixplanet_tools

## take 10000 random seuqnces for testing
##module load seqkit/2.2.0
##seqkit sample -n 10000 /g/scb2/bork/ustick/nfixplanet_pipeline/00_data_files/sequence_files/nfix_gene_seqs.fna > /g/scb2/bork/ustick/nfixplanet_pipeline/00_data_files/sequence_files/nfix_gene_seqs_10000.fna

## 01 ## run prodigal on genomes then merge together reads ##
#./01.01.00_prodigal.sh

## 02 ## run hmmer on sequences ##
#sbatch 01.02.00_run_hmms.sh /g/scb2/bork/ustick/nfixplanet_pipeline/00_data_files/genomes/prodigal_calls/01_all_genes.fna
out_path="/g/scb2/bork/ustick/nfixplanet_pipeline/01_annotation_pipeline/data/ustick_test_output"
## update to simplified format
#./01.02.01_format_hmmer_tbl.py $out_path/02_test_hmm.tbl > $out_path/02_custom_results.tsv
## take only top hit ##
## deactivate conda for the rscript ##
#./01.02.02_only_best_hit.R $out_path/02_custom_results.tsv $out_path/02_top_hit.tsv

## 03 ## qc based on threshold & make sure they are all on same contig ##
#./01.03.00_extract_nfix.sh

## 04 ## qc based on genomic context ##
#./01.04.00_neighborhood_nifanfvnf.sh nif 10
#./01.04.00_neighborhood_nifanfvnf.sh anf 10
#./01.04.00_neighborhood_nifanfvnf.sh vnf 10
#./01.04.00_neighborhood_Chl.sh 10
#./01.04.00_neighborhood_nfl.sh 10

