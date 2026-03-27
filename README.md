# Snakemake-based Variant Annotation Pipeline

## Overview
This pipeline is a Snakemake-based workflow that generates Mutation Annotation Format (MAF) files from VCF files and performs variant classification.

## Features
- Automatically converts VCF files to MAF format
- Integrates multiple annotation tools
- Performs variant classification and VUS prioritization
- Ensures reproducibility using Snakemake

## Dependencies
Please install the following tools and data manually:

- InterVar  
- ANNOVAR  
- oncokb-annotator  
- Apptainer (Singularity)
- VCF file of ClinVar
- VEP cache data
- VCF file of Population allele frequency (e.g. jMorp)

## Environment
- OS: Rocky Linux v9.6  
- Snakemake: v9.13.7  

## Installation

```bash
git clone https://github.com/hnakahara/S-VAP.git
cd S-VAP
```

## Setup

### 1. Edit config.yaml

Modify the following files in the `config/` directory according to your environment:

- `config.yaml` (GRCh37)
- `config_hg38.yaml` (GRCh38)

Parameters:

- reference: Path to the reference FASTA file  
- cachedir: Path to VEP cache directory (default: v112)  
- tommo: Population frequency file  
- OncoDir: Directory containing oncokb-annotator  
- token: OncoKB API token  
- intervar: Directory where InterVar is installed  
- annovar: Directory where ANNOVAR is installed  
- clinvar: Path to ClinVar VCF file  
- spliceai_snv: Precomputed SpliceAI SNV file  
- spliceai_indel: Precomputed SpliceAI INDEL file  

### 2. Edit vcfs.tsv

Edit `config/vcfs.tsv`:

| Column | Description |
|--------|-------------|
| samples | Sample identifier |
| normal | Path to input VCF file |

## Usage

### Basic execution

```bash
snakemake --use-conda --use-apptainer --cores 4 --apptainer-args --rerun-triggers mtime
```

### With bind option

```bash
snakemake --use-conda --use-apptainer --cores 4 --apptainer-args "--bind /your_directory" --rerun-triggers mtime
```

### Troubleshooting

If a VEP plugin error occurs during the `vcf2maf` step, you may need to replace the `vcf2maf.pl` script.

You can download the latest version from the official repository:
https://github.com/mskcc/vcf2maf

## Output

A `results/` directory will be generated under `S-VAP/`, containing the following subdirectories:

- 08vcf2maf
- 09oncokb-annotator
- 10intervar
- 11final-output

Each directory contains per-sample annotation results.

### Final output (11final-output)

- Per-sample merged results
- Combined MAF file for all samples

#### Output files

- merge_raw.maf  
  → Before variant classification  

- merge_variant_class.maf  
  → After variant classification  

- merge_variant_class_vus_prior.maf  
  → After VUS prioritization  

## MAF Format

Refer to:
https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

## Optional: Generating VCF files

If VCF files are not available, they can be generated using TransVar.

### Environment setup

```bash
conda create -n transvar python=3.9 conda-forge::pandas -y
conda activate transvar
pip install transvar
```

### Configuration

```bash
export TRANSVAR_CFG=/path/to/transvar.cfg
```

### Run

```bash
python CreateVcf.py \
  --input /path/to/input.csv \
  --outdir /path/to/output \
  --prefix sample
```

**Input example:**  
Refer to `example/input_for_transvar.csv` for the expected input format.

Example format (tab-delimited):

| gene  | cdna      | protein | transcript     |
|-------|-----------|--------|----------------|
| MUTYH | 892-2A>G  |        | NM_001048171   |
| ABL1  | 925C>G    | P309A  | NM_005157      |
| TERT  | -124C>T   |        | NM_198253      |

- The header names must not be changed.
- The `protein` column is optional. It can be filled or left blank and does not affect the processing.


## License

This project is licensed under the Creative Commons Attribution-NonCommercial (CC-BY-NC) license.  
For more details, refer to the LICENSE file.

## Web Application (Beta) and Public API Demo

A beta version of the web application implementing a portion of the current pipeline algorithm is available to a limited set of users in Japan.

For demonstration purposes, a publicly accessible REST API is also provided:

- API documentation (Swagger UI):  
  https://1603-027.a.hiroshima-u.ac.jp/api/swagger/