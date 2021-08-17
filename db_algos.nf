#!/usr/bin/env nextflow


process sortmerRNA_SE{
	    
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    scratch params.scratch_small
    publishDir path: "${params.WD}/SortmeRNA/"    

    input:
	path trimmed_fastq_SE
        val SortmeRNA_Reflist
    
    output:
	path "${trimmed_fastq_SE.baseName}_aligned*"
	path "${trimmed_fastq_SE.baseName}_aligned.log", emit: sortmerna_log
    	path "${SE_mRNA_read}.fastq", emit: SE_mRNA_read
	
    script:
	SE_mRNA_read =  trimmed_fastq_SE.baseName.replace('trim_','')   
        SortmeRNA_Reflist = SortmeRNA_Reflist.collect{'--ref '+ it }.join(' ')

    
"""

    mkdir kvdb

    sortmerna \
        ${SortmeRNA_Reflist} \
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




process bowtie2_build{

    cpus params.htp_cores 
    memory "${params.m_mem} GB"    
    publishDir path: "${bt2_index_path}"
    input:
	path fasta_references
    	val bt2_index_base
    	val bt2_index_path

    output:
	path "${bt2_index_base}*"


"""

    bowtie2-build \
	${fasta_references} \
        ${bt2_index_base}
                        		      
"""

}




process bowtie2_SE{

    cpus params.htp_cores 
    memory "${params.m_mem} GB"
    publishDir path: "${params.WD}/bowtie2"
    input:
        path SE_reads
	val fasta_references
	val bt2_index_base
    	val bt2_index_path

    output:
	 path SAM_file, emit: SAM
	 path metrics
	 path quick_stats

script:
SAM_file =  SE_reads.getSimpleName() + '.sam'
metrics = SE_reads.getSimpleName() + '.metrics'
quick_stats = SE_reads.getSimpleName() + '.quick_stats'

"""
  
   bowtie2 \
      --threads ${params.htp_cores} \
      -q \
      --no-unal \
      -k 20 \
      -x ${bt2_index_path}/${bt2_index_base} \
      -U ${SE_reads} \
      -S ${SAM_file} \
      --met-file ${metrics} 2> ${quick_stats} 
      
"""	

}








