#!/usr/bin/env nextflow




process sortmerRNA_SE{
	    
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    scratch params.scratch_small
    publishDir path: "${params.WD}/SortmeRNA/"
    
    input:
	path trimmed_fastq_SE
	val SortmeRNA_idx_dir 
        val SortmeRNA_Reflist
    
    output:
	path "${trimmed_fastq_SE.baseName}_aligned*"
	path "${trimmed_fastq_SE.baseName}_aligned.log", emit: sortmerna_log
    	path "${SE_mRNA_read}.fastq", emit: SE_mRNA_read
	
    script:
	SE_mRNA_read =  trimmed_fastq_SE.baseName.replace('trim_','')   
        SortmeRNA_Reflist = SortmeRNA_Reflist.collect{'--ref '+ it }.join(' ')

    
"""
    
    mkdir -p "${HOME}/sortmerna/run/"
    mkdir kvdb
    
    sortmerna \
        ${SortmeRNA_Reflist} \
    	--reads ${trimmed_fastq_SE} \
    	--aligned ${trimmed_fastq_SE.baseName}_aligned \
        --kvdb kvdb \
        --idx-dir  ${SortmeRNA_idx_dir} \
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



process makeblastdb{

    cpus params.htp_cores 
    memory "${params.m_mem} GB"    
    publishDir path: "${DB_path}"
    input:
	path fasta_reference
    	val molecule_type
	val DB_path

    output:
	path "${fasta_reference}*"

"""

    makeblastdb \
        -in ${fasta_reference} \
        -dbtype ${molecule_type}
                        		      
"""


} 



process blast_tophit{

    cpus params.htp_cores 
    memory "${params.m_mem} GB"    
    publishDir path: "${params.WD}/Assembly"
    input:
	path query_seqs
	path fasta_reference
        val DB_path
    output:
	path blastout_fmt6
    
script:
db_name =  fasta_reference.getBaseName()   	
blastout_fmt6 = query_seqs.getSimpleName() + '.blastx'
"""

    blastx \
       -query ${query_seqs} \
       -db  ${DB_path}/${fasta_reference} \
       -out ${blastout_fmt6} \
       -evalue 1e-20 \
       -num_threads ${params.htp_cores} \
       -max_target_seqs 1\
       -outfmt 6 
                   		      
"""


} 



process busco_auto_euk{

    cpus params.mtp_cores 
    memory "${params.m_mem} GB"
    publishDir path: "${params.WD}"
    input:
        path denovo_ref
    	val species

    output:
        path "Busco_${denovo_ref}"


"""
  
    busco \
        -m transcriptome \
        -i ${denovo_ref} \
        --species ${species} \
        --cpu ${params.mtp_cores} \
        -o Busco_${denovo_ref} \
        --download_path ${DB_REF} \
        --auto-lineage-euk
        

"""	

}





