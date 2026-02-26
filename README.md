# NFixPlanet

Python package for detection and quantification of nitrogen-fixing microorganisms (diazotrophs) from genomes and short-read metagenomes.

## Description

NFixPlanet provides workflows for identifying nitrogen fixation genes in genomes/contigs and quantifying the abundance and taxonomic composition of diazotrophs in metagenomic datasets. It combines profile Hidden Markov Model (HMM) annotation, genomic context validation, and coverage-based abundance estimation using a curated diazotroph reference database.

The package contains two main workflows:
- **Genome annotation** (`annotate`) — identifies nitrogen fixation genes and operons in genome assemblies using HMMs and genomic context filtering.
- **Metagenome quantification** (`profile`) — maps short reads to a diazotroph reference database and computes gene abundance, genome abundance, and taxonomic relative abundance.

The metagenome workflow currently supports short-read sequencing data only, as abundance estimation relies on short-read coverage profiling.

## Installation
The recommended method of installation is via [bioconda](https://anaconda.org/channels/bioconda/packages/nfixplanet/overview)
```bash
conda install -c bioconda nfixplanet
```

### Manual installation
Alternatively, Nfixplanet can be installed via `pip` and the requiements can be installed manually
```bash
pip install nfixplanet
```

Requirements:
- Python 3.11
- `defopt < 7`
- `wget`
- [Prodigal v2.6.3](https://github.com/hyattpd/Prodigal)
- [HMMER v3.4](http://hmmer.org/)
- [fastp v0.24.0](https://github.com/OpenGene/fastp)
- [hostile v2.0.0](https://github.com/bede/hostile)
- [minimap2 v2.28](https://github.com/lh3/minimap2)
- [coverm v0.7.0](https://github.com/wwood/CoverM)


### Local installation

```bash
git clone git@git.embl.org:grp-bork/nfixplanet.git
cd nfixplanet
conda create -c bioconda -n nfixdev python=3.11 prodigal=2.6.3 hmmer=3.4 "defopt<7" wget hostile=2.0.0 fastp=0.24.0 minimap2=2.28 coverm=0.7.0
conda activate nfixdev
pip install -e .[dev]
```

## Usage
### nfixplanet annotate

Pipeline for identifying nitrogen fixation genes and operons in genome assemblies.

This workflow:
1. Predicts open reading frames (ORFs) using Prodigal (optional if ORFs are provided)
2. Searches ORFs against curated nitrogen fixation HMM profiles using HMMER
3. Retains the best HMM hit per ORF
4. Filters results to ensure required genes occur on the same contig
5. Performs genomic context validation based on operon structure and gene proximity
6. Resolves ambiguous gene assignments (e.g., nifH vs vnfH) using neighboring genes

Output consists of high-confidence nitrogen fixation gene annotations and operon assignments.

This workflow is designed for any nucleotide sequence including genome assemblies, metagenome-assembled genomes (MAGs), or individual contigs.

Basic command:

```bash
nfixplanet annotate --input_fasta /path/to/fasta --output_directory /path/to/output
```
Optional arguments:

- `--genomic_context_range <int>`: Maximum number of genes upstream or downstream to consider for operon context (default: 10).
- `--cpus <int>`: Number of CPUs to use for HMMscan (default: 2, max recommended: 4).
- `--verbose`: Enable verbose (DEBUG) logging.
- `--version`: Print version number and exit.

### nfixplanet profile
Pipeline for quantifying diazotroph gene, genome, and taxonomic abundance from short-read metagenomes.

This workflow integrates read mapping, abundance estimation, and taxonomic profiling into a single pipeline.

Steps include:
1. Pre-processing quality control on metagenomic reads
2. Mapping metagenomic reads to a curated diazotroph reference gene database using CoverM
3. Calculating coverage per gene across samples
4. Aggregating gene coverage into genome-level abundance estimates
5. Normalizing genome abundance by gene count per genome
6. Assigning genomes to taxonomy using a reference taxonomy table
7. Aggregating genome abundance into taxonomic relative abundance profiles

Outputs include:
1. Gene abundance table (coverage per nitrogen fixation gene per sample)
2. Genome abundance table (normalized coverage per diazotroph genome per sample)
3. Taxonomic relative abundance tables at multiple taxonomic levels (e.g., phylum, class, order, family, genus)

These outputs provide both functional and taxonomic quantification of diazotroph communities across metagenomic samples.

Basic command:

```bash
nfixplanet profile \
  --sample_id SAMPLE_NAME \
  --read_1 /path/to/read_1.fastq \
  --read_2 /path/to/read_2.fastq \
  --single /path/to/reads.fastq \
  --output_directory /path/to/output
```

Required arguments:
- `--sample_id <str>`: Name of the FASTA/FASTQ sample.
- `--output_directory <str>`: Path to output directory.

Input read options (choose one mode):
- `--read_1 <str>`: Path to FASTA/FASTQ file for paired-end read 1 (R1).
- `--read_2 <str>`: Path to FASTA/FASTQ file for paired-end read 2 (R2). Must be provided together with `--read_1`.
- `--single <str>`: Path to FASTA/FASTQ file for single-end reads Can be used on its own or with `--read_1` and `--read_2`.

Optional arguments:
- `--work_directory <str>`: Path to directory for temporary files (default: tmp).
- `--cpus <int>`: Number of CPUs used by processes (default: 8).
- `--verbose`: Enable verbose logging.

## Authors and acknowledgment
- [Mahdi Robbani](https://github.com/mahdi-robbani)
- [Lucas Ustick](https://github.com/ljustick)
- [Anthony Fullam](https://github.com/fullama)

## License
This project is licensed under the MIT License. See the LICENSE file for details.
