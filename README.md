# BIOF501A Term Project: Taxonomic Classification and Visualization of Short-Read Metagenomic Sequencing Data

## Background and Rationale

## Usage
```
# download repository
git clone https://github.com/AmosFong1/BIOF501A

# create conda environment from environment.yml
conda env create -f environment.yml

# download KrakenTools
git clone https://github.com/jenniferlu717/KrakenTools

# remove faulty sym link from krona install
rm -rf "$(pwd)"/.conda/envs/BIOF501A/opt/krona/taxonomy

# make directory to store krona database
mkdir -p "$(pwd)"/krona/taxonomy

# create sym link to krona database
ln -s "$(pwd)"/krona/taxonomy "$(pwd)"/.conda/envs/BIOF501A/opt/krona/taxonomy

# download krona database
wget -pO "$(pwd)"/krona/taxonomy/taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# run ktUpdateTaxonomy script
ktUpdateTaxonomy.sh --only-build

# run nextflow script
nextflow run BIOF501A.nf
```

## Input
The workflow inputs include raw metagenomic sequencing data (FASTQ) from a global sewage-based antimicrobial resistance (AMR) profiling project. The original analysis was published in a paper titled "Genomic analysis of sewage from 101 countries reveals global landscape of antimicrobial resistance". The FASTQ files generated from this study can be accessed at the European Nucleotide Archive under accession numbers PRJEB40798, PRJEB40816, PRJEB40815, PRJEB27621, PRJEB51229, and ERP015409. The sample used to model this workflow was collected on 2017-06-23 in Vancouver, BC. The FASTQ files for this sample can be found under study accession PRJEB27621, sample accession SAMEA4777410, experiment accession ERX2697767, and run accession ERR2683153. The other inputs include the Bowtie2 database of the human host genome GRCh38 (hg38), and the Kraken 2 / Bracken RefSeq indexes (Standard-8 collection), which are both downloaded and unpackaged as part of the workflow.

## Output
The workflow's main output is the `krona.html` file, which can be found in the `data/krona` directory. This file provides an interactive metagenomic visualization of estimated taxonomical abundances that can be downloaded and explored with any web browser. The other outputs include the `fastp_ERR2683153.fastp.html`, which can be found in the `data/fastp` directory, and the `fastp_ERR2683153_1_fastqc.html` and `fastp_ERR2683153_2_fastqc.html` files, which can be found in the `data/fastqc` directory. Together, these files provide a report of the FASTP pre-processing run and the FASTP-processed FASTQ files, which can be downloaded and explored with any web browser.
