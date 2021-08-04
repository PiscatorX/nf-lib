#!/usr/bin/env nextflow

process sortmerRNA_SE{

    
    //cpus params.ltp_cores
    publishDir path: "${params.WD}/SortmeRNA"
    memory "${params.m_mem} GB"

    input:
	path trimmed_fastq_SE
        val SILVA
    
    output:
	path "${trimmed_fastq_SE.baseName}_aligned*"
	path "${trimmed_fastq_SE.baseName}_aligned.log", emit: sortmerna_log
    	path "${SE_mRNA_read}.fastq", emit: SE_mRNA_read
	
script:
	SE_mRNA_read =  trimmed_fastq_SE.baseName.replace('trim_','')   


"""

    mkdir kvdb

    sortmerna \
        --ref ${params.SILVA} \
	--reads ${trimmed_fastq_SE} \
	--aligned ${trimmed_fastq_SE.baseName}_aligned \
        --kvdb kvdb \
	--other mRNA \
	--threads ${params.htp_cores} \
	--fastx \
	-m ${params.h_mem}000 \
	--task 4 \
	-v 

      mv mRNA.fq ${SE_mRNA_read}.fastq   
   
"""
//underscores are sometimes removed in read name
}



