# Taxonomic Classification and Visualization of Short-Read Metagenomic Sequencing Data
### By: Amos Fong

***

## Background and Rationale

Global genomic pathogen surveillance stands at the forefront of public health initiatives, providing an essential tool that enables public health scientists to identify and address the emergence of pathogenic strains before they escalate from an endemic stage to a global pandemic. The monitoring of urban sewage represents a pragmatic and ethically sound approach, serving as a surrogate for tracking the epidemiological dynamics and selection of both bacterial strains and their antimicrobial resistance genes. Urban sewage, as a sampling medium, provides a comprehensive snapshot of the microbiological contributions from both local humans and animals, as well as contributions from the immediate environment. Notably, it eliminates the need to request consent for collecting human fecal samples and bypasses the laborious effort of assembling a cohort large enough to generalize findings to the population level.

Shotgun metagenomics sequencing offers a high-throughput and hypothesis-naive approach to capturing the complex genetic diversity inherent in environmental samples, such as urban sewage. A key step in the analysis of metagenomic sequencing data is the taxonomic profiling of sequencing reads to infer the sample's microbial makeup. Among the currently available bioinformatic tools used for taxonomic classification, Kraken and its derivative tools Braken/KrakenUniq/Kraken2 stand out as some of the best-performing tools in terms of computational requirements and performance accuracy. Furthermore, the inherent uncertainty in taxonomic classifications underscores the importance of considering their hierarchical contexts and prediction confidence scores when visualizing metagenomic data. The Krona package represents an attractive solution for metagenomics data visualization, capturing both the quantitative hierarchical relationships between bacterial species and enabling the visualization of secondary metadata variables.

Here, we present a bioinformatics workflow that performs end-to-end quality control, read pre-processing, host decontamination, taxonomic classification, abundance estimation, and interactive visualization for short-read metagenomics sequencing data. We anticipate that this workflow will provide a comprehensive solution for querying and analyzing metagenomics sequencing data from publicly available FASTQ files deposited in the NCBI Sequence Read Archive (SRA). The design of this workflow was inspired by the [nf-core taxaprofiler]( https://github.com/nf-core/taxprofiler) pipeline.

This workflow begins by downloading and unpacking the Bowtie2 database of the human host genome GRCh38 (hg38) and the Kraken 2 / Bracken RefSeq indexes (Standard-8 collection). Next, the workflow loads in raw metagenomic sequencing data in FASTQ format from the NCBI SRA. For the purposes of modelling this workflow, we will be analyzing sequencing results from a urban sewage sample collected in Vancouver, BC. Subsequently, the workflow utilizes FASTP to perform read pre-processing. FASTP is a tool designed to provide a fast all-in-one preprocessing step for FASTQ files. The workflow then uses FASTQC to generate read quality reports, providing a summary of the FASTQ read quality after FASTP pre-processing. In parallel, Bowtie2 is used to remove reads mapping to the hg38 genome found in the FASTP pre-processed reads. Following these steps, the workflow utilizes Kraken2 to assign taxonomic labels to the Bowtie2-filtered, FASTP pre-processed reads. Subsequently, Bracken is used to estimate the taxonomic abundances at the species level using the Kraken2-assigned classes. The workflow then uses the the python script `kreport2krona.py` from KrakenTools to convert the Kraken2/Bracken report into a Krona-compatible text file. Finally, the workflow uses Krona to generate a HTML file which provides an interactive visualization of the estimated taxonomical abundances.

![workflow](dag.svg)

## SOP
### Dependencies
To run this workflow, the user must have [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed. Additionally, this workflow depends on the following packages:
```
bowtie2=2.5.2
bracken=2.9
fastp=0.23.4
fastqc=0.12.1
kraken2=2.1.3
krakentools=1.2
krona=2.8.1
nextflow=23.10.0
```
### Installation
Step 1: Deactivate conda environment
```
conda deactivate
```
Step 2: Clone repository
```
git clone https://github.com/AmosFong1/BIOF501A
```
Step 3: Navigate to project directory
```
cd BIOF501A
```
Step 4: Create conda environment
```
conda env create -f environment.yml
```
Step 5: Clone KrakenTools repository
```
git clone https://github.com/jenniferlu717/KrakenTools
```
Step 6: Remove faulty sym link
```
rm -rf "$(pwd)"/.conda/envs/BIOF501A/opt/krona/taxonomy
```
Step 7: Create directory to store new krona database
```
mkdir -p "$(pwd)"/krona/taxonomy
```
Step 8: Create sym link to new krona database
```
ln -s "$(pwd)"/krona/taxonomy "$(pwd)"/.conda/envs/BIOF501A/opt/krona/taxonomy
```
Step 9: Download new krona database
```
wget -pO "$(pwd)"/krona/taxonomy/taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
```
Step 10: Run `ktUpdateTaxonomy.sh` script
```
ktUpdateTaxonomy.sh --only-build
```
### Usage
Step 1: Activate conda environment
```
conda activate BIOF501A
```
Step 2: Run the workflow (use `-resume` option as needed)
```
nextflow run BIOF501A.nf
nextflow run BIOF501A.nf -resume
```

## Input
The workflow inputs include raw metagenomic sequencing data (FASTQ) from a global sewage-based antimicrobial resistance (AMR) profiling project. The original analysis was published in a paper titled "Genomic analysis of sewage from 101 countries reveals global landscape of antimicrobial resistance". The FASTQ files generated from this study can be accessed at the European Nucleotide Archive under accession numbers PRJEB40798, PRJEB40816, PRJEB40815, PRJEB27621, PRJEB51229, and ERP015409. The sample used to model this workflow was collected on 2017-06-23 in Vancouver, BC. The FASTQ files for this sample can be found under study accession PRJEB27621, sample accession SAMEA4777410, experiment accession ERX2697767, and run accession ERR2683153. The other inputs include the Bowtie2 database of the human host genome GRCh38 (hg38), and the Kraken 2 / Bracken RefSeq indexes (Standard-8 collection), which are both downloaded and unpackaged as part of the workflow.

## Output
The workflow's main output is the `krona.html` file, which can be found in the `data/krona` directory. This file provides an interactive metagenomic visualization of estimated taxonomical abundances that can be downloaded and explored with any web browser. The other outputs include the the `fastp_ERR2683153_1_fastqc.html` and `fastp_ERR2683153_2_fastqc.html` files, which can be found in the `data/fastqc` directory. These files provide a QC report of the FASTP-processed FASTQ files, which can be downloaded and explored with any web browser.

![output](krona.svg)
