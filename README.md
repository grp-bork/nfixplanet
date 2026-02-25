# NFixPlanet

Python package for detection of nitrogen fixers.

## Description
TODO: need to update description
TODO: mention somewhere that the pipeline only works for short read sequences (coverm step)

## Installation
TODO: install via conda
The recommended method of installation is via [bioconda](TODO: add to conda package)
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

Pipeline for annotating genomes as N-fixers TODO: update description

Basic command:

```bash
nfixplanet annotate --input_fasta /path/to/fasta --output_directory /path/to/output
```
Optional arguments:

- `--genomic_context_range <int>`: Maximum number of genes upstream or downstream to consider for operon context (default: 10).
- `--cpus <int>`: Number of CPUs to use for HMMscan (default: 2, max recommended: 4).
- `--verbose`: Enable verbose (DEBUG) logging.
- `--version`: Print version number and exit.

### nfixplanet map
Pipeline for mapping metagenomes TODO: update description

Basic command:

```bash
nfixplanet map \
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
- Lucas Ustick
- [Anthony Fullam](https://github.com/fullama)

## License
This project is licensed under the MIT License. See the LICENSE file for details.
