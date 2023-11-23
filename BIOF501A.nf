#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastp = "${baseDir}/data/fastq/*_{1,2}.fastq"
params.fastqc = "${baseDir}/data/fastp/fastp_*_{1,2}.fastq"
params.bowtie2 = "${baseDir}/data/fastp/fastp_*_{1,2}.fastq"
params.kraken2 = "${baseDir}/data/bowtie2/fastp_*_host_removed_{1,2}.fastq.gz"

workflow {

DOWNLOAD_READS()

DOWNLOAD_HOST()

DOWNLOAD_DATABASE()

fastp_ch = Channel.fromFilePairs(params.fastp)
FASTP(fastp_ch)

fastqc_ch = Channel.fromPath(params.fastqc)
FASTQC(fastqc_ch)

bowtie2_ch = Channel.fromFilePairs(params.bowtie2)
BOWTIE2(bowtie2_ch)

kraken2_ch = Channel.fromFilePairs(params.kraken2)
bracken_ch = KRAKEN2(kraken2_ch)
kreport2krona_ch = BRACKEN(bracken_ch)
krona_ch = KREPORT2KRONA(kreport2krona_ch)
KRONA(krona_ch)

}

process DOWNLOAD_READS {
    publishDir = "${baseDir}/data/fastq"
    output:
        path("ERR2683153_{1,2}.fastq")
    script:
        """
        wget -qO ERR2683153_1.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/003/ERR2683153/ERR2683153_1.fastq.gz && gunzip -c ERR2683153_1.fastq.gz > ERR2683153_1.fastq
        wget -qO ERR2683153_2.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/003/ERR2683153/ERR2683153_2.fastq.gz && gunzip -c ERR2683153_2.fastq.gz > ERR2683153_2.fastq
        """
}

process DOWNLOAD_HOST {
    publishDir "${baseDir}/data/host"
    output:
        path("GRCh38_noalt_as")
    script:
        """
        wget -q https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip && unzip GRCh38_noalt_as.zip
        """
}

process DOWNLOAD_DATABASE {
    publishDir "${baseDir}/data/database"
    output:
        path("hash.k2d")
        path("opts.k2d")
        path("taxo.k2d")
        path("seqid2taxid.map")
        path("inspect.txt")
        path("database100mers.kmer_distrib")
        path("database150mers.kmer_distrib")
        path("database200mers.kmer_distrib")
        path("database250mers.kmer_distrib")
        path("database300mers.kmer_distrib")
        path("database50mers.kmer_distrib")
        path("database75mers.kmer_distrib")
        path("ktaxonomy.tsv")
        path("library_report.tsv")
        path("unmapped_accessions.txt")
    script:
        """
        wget -q https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz && tar -xvzf k2_standard_08gb_20231009.tar.gz
        """
}

process FASTP {
publishDir "${baseDir}/data/fastp"
input:
    tuple val(sample), path(reads)
output:
    path("fastp_${sample}_{1,2}.fastq")
    path("fastp_${sample}.fastp.html")
script:
def (r1, r2) = reads
    """
    fastp \\
        --in1 "${r1}" \\
        --in2 "${r2}" \\
        --out1 "fastp_${sample}_1.fastq" \\
        --out2 "fastp_${sample}_2.fastq" \\
        --html "fastp_${sample}.fastp.html"
    """
}

process FASTQC {
    publishDir "${baseDir}/data/fastqc"
    input:
        path(reads)
    output:
        path("*.html")
    script:
    """
    fastqc ${reads}
    """
}

process BOWTIE2 {
    publishDir "${baseDir}/data/bowtie2"
    input:
        tuple val(sample), path(reads)
    output:
        path("*_host_removed_{1,2}.fastq.gz")
    script:
    def (r1, r2) = reads
    """
    idx_base=\$(find ${baseDir}/data/host/GRCh38_noalt_as/ -name '*.bt2' | awk -F \".\" '{print \$1 | \"sort -u\"}')
    
    bowtie2 -p 8 -x \${idx_base} \
        -1 "${r1}" \
        -2 "${r2}" \
        --un-conc \
        ${sample}_host_removed \
        > ${sample}_mapped_and_unmapped.sam

    mv ${sample}_host_removed.1 ${sample}_host_removed_1.fastq.gz
    mv ${sample}_host_removed.2 ${sample}_host_removed_2.fastq.gz
    """
}

process KRAKEN2 {
    publishDir "${baseDir}/data/kraken2"
    input:
        tuple val(sample), path(reads)
    output:
        path("k2_report.txt")
    script:
    def (r1, r2) = reads
    """
    kraken2 --db ${baseDir}/data/database \\
        --report k2_report.txt \\
        --output k2_output.txt \\
        --paired "${r1}" "${r2}"
    """
}

process BRACKEN {
    publishDir "${baseDir}/data/bracken"
    input:
        path(report)
    output:
        path("b_report.txt")
    script:
    """
    bracken -d ${baseDir}/data/database -i ${report} -r 100 -l S -t 10 -o b_output.txt -w b_report.txt
    """
}

process KREPORT2KRONA {
    publishDir "${baseDir}/data/kreport2krona"
    input:
        path(report)
    output:
        path("krona.txt")
    """
    python ${baseDir}/KrakenTools/kreport2krona.py -r ${report} -o krona.txt --no-intermediate-ranks
    """
}

process KRONA {
    publishDir "${baseDir}/data/krona"
    input:
        path(report)
    output:
        path("krona.html")
    """
    ktImportText ${report} -o krona.html
    """
}