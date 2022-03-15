[![Continuous Integration](https://github.com/LUMC/haploblock-combiner/actions/workflows/ci.yml/badge.svg)](https://github.com/LUMC/haploblock-combiner/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
![GitHub release](https://img.shields.io/github/v/release/LUMC/haploblock-combiner)
![Commits since latest release](https://img.shields.io/github/commits-since/LUMC/haploblock-combiner/latest)

# haploblock-combiner
Generate all possible combinations between phased haploblocks.

## Installation
Download the repository from github
```bash
git clone https://github.com/LUMC/haploblock-combiner.git
```

Install and activate the
[conda](https://docs.conda.io/en/latest/miniconda.html)
environment.
```bash
conda env create --file environment.yml
conda activate haploblock-combiner
```

## Settings
There are three levels where configuration options are set, in decreasing order
of priority.
1. Flags passed to snakemake using `--config`, or in the specified
   `--configfile`.
2. Setting specified in the PEP project configuration, under the key
   `haploblock-combiner`
3. The default settings for the pipeline, as specified in the `common.smk` file

## Usage
The following settings are required:
`--config pepfile=pepfile.yml reference=/path/to/ref.fa region=chrM`

**Note:** Please make sure to specify as small a region as possible. This
pipeline will generate 2^(n-1) files, where n is determined by the number of
phase sets in the input VCF. So for an unphased VCF with 20 heterozygous
variants, it will generate over **1 million files**.

## Limitations
* Only VCF files phased by
[WhatsHap](https://whatshap.readthedocs.io/en/latest/) are currently supported.
* Only VCF files that are compressed with bgzip and indexe with tabix are
currently supported
