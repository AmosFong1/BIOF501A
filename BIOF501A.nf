#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

DOWNLOAD_HOST()
DOWNLOAD_DATABASE()

ids = ["ERR2683153"]
ch_reads = Channel.fromSRA(ids)
FASTP(ch_reads)
ch_fastp = FASTP.out.reads
FASTQC(ch_fastp)
BOWTIE2(ch_fastp)
ch_bowtie2 = BOWTIE2.out.reads
ch_bracken = KRAKEN2(ch_bowtie2)
ch_kreport2krona = BRACKEN(ch_bracken)
ch_krona = KREPORT2KRONA(ch_kreport2krona)
KRONA(ch_krona)

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
    tuple val(sample), path("fastp_${sample}_{1,2}.fastq"), emit: reads
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
        tuple val(sample), path(reads)
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
        tuple val(sample), path("*_host_removed_{1,2}.fastq.gz"), emit: reads
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