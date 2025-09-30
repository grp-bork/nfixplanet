# NFixPlanet

Python package for detection of nitrogen fixers.

## Description
TODO: need up update description

## Installation
TODO: install via pip
TODO: install via conda

### Local installation

```
git clone git@git.embl.org:grp-bork/nfixplanet.git
cd nfixplanet
conda create -n nfixtest python=3.10
conda activate nfixtest
conda install prodigal hmmer
pip install -e .
```

Requirements:
- Prodigal V2.6.3: February, 2016
- HMMER 3.4 (Aug 2023); http://hmmer.org/

## Usage

Basic command:

```bash
nfixplanet --input_fasta /path/to/fasta --output_directory /path/to/output
```
Optional arguments:

- `--genomic_context_range <int>`: Maximum number of genes upstream or downstream to consider for operon context (default: 10).
- `--cpus <int>`: Number of CPUs to use for HMMscan (default: 2, max recommended: 4).
- `--verbose`: Enable verbose (DEBUG) logging.
- `--version`: Print version number and exit.

## Authors and acknowledgment
- [Mahdi Robbani](https://github.com/mahdi-robbani)
- Lucas Ustick
- [Anthony Fullam](https://github.com/fullama)

## License
This project is licensed under the MIT License. See the LICENSE file for details.
