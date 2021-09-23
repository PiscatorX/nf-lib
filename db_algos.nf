#!/usr/bin/env nextflow




process sortmerRNA_SE{

    publishDir path: "${params.WD}/SortmeRNA/"
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    scratch params.scratch_small
    input:
	path trimmed_fastq_SE
	val SortmeRNA_idx_dir 
        val SortmeRNA_Reflist
    
    output:
	path "${trimmed_fastq_SE.baseName}_aligned.fastq"
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
    	--aligned ${trimmed_fastq_SE.baseName}_aligned.fastq \
        --kvdb kvdb \
        --idx-dir  ${SortmeRNA_idx_dir} \
    	--other mRNA \
    	-a ${params.ltp_cores} \
    	--fastx \
    	-m ${params.l_mem}000 \
    	--task ${params.mtp_cores} \
    	-v 

      mv mRNA.fq ${SE_mRNA_read}.fastq   
       
"""
//underscores are sometimes removed in read name
}




process bowtie2_build{

    storeDir "${bt2_index_path}"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.mtp_cores 
    input:
	path fasta_reference
    	val bt2_index_path

    output:
	path "${index_basename}*"

script:
index_basename = fasta_reference.getName()
"""

    bowtie2-build \
        --threads ${params.mtp_cores} \
	${fasta_reference} \
        ${index_basename}
                        		      
"""

}



process bowtie2_SE{

    publishDir path: "${params.WD}/bam/", pattern: "*.bam"
    publishDir path: "${params.WD}/sam/", pattern: "*.sam"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.htp_cores     
    
    
    input:
        path SE_reads
	path fasta_reference
    	val bt2_index_path

    output:
	path SAM_file, emit: SAM
	path BAM_file, emit: BAM
	path metrics
	path quick_stats

script:
SAM_file =  SE_reads.getSimpleName() + '.sam'
BAM_file =  SE_reads.getSimpleName() + '.bam'
metrics = SE_reads.getSimpleName() + '.metrics'
quick_stats = SE_reads.getSimpleName() + '.quick_stats'
index_basename = fasta_reference.getName()

"""
  
   bowtie2 \
      --threads ${params.mtp_cores} \
      -x ${bt2_index_path}/${index_basename} \
      -U ${SE_reads} \
      -S ${SAM_file} \
      --met-file ${metrics} 2> ${quick_stats} 
         
   
   samtools \
       view \
       ${SAM_file} \
       -F 4 \
       -b \
       --threads ${params.mtp_cores} \
       -o ${BAM_file} 
            
"""	

}




process bam_index{

    publishDir path: "${params.WD}/indexed_bam/"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.mtp_cores 
       
    input:
        path BAM
	val fasta_reference

    output:
        path "${BAM.baseName}*", emit: sorted_BAM

"""

    samtools \
        sort \
        ${BAM} \
        -m "${params.l_mem}G" \
        --threads  ${params.mtp_cores} \
        --reference ${fasta_reference} \
        -o  ${BAM.baseName}_sorted.bam

    samtools \
        index \
        ${BAM.baseName}_sorted.bam \
        -@  ${params.mtp_cores} 


"""
      
}




process samtools_fasta_ref{

    publishDir "${params.WD}/indexed_bam/"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.htp_cores

    input:
        path fasta_reference
	val  name

     output:
        path "*.fai"
	path "${name}.fasta"

script:
fname = fasta_reference.getName()
"""
   
    mv ${fname}  ${name}.fasta 	

    samtools \
       faidx \
       ${name}.fasta 
       
      
"""

}




process makeblastdb{

    publishDir path: "${DB_path}"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.htp_cores
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

    publishDir "${params.WD}/Blast_top_hits/"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.mtp_cores
    
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
       -num_threads ${params.mtp_cores} \
       -max_target_seqs 1\
       -outfmt 6 
                   		      
"""


} 



process busco_auto_euk{

    publishDir path: "${params.WD}"
    scratch params.scratch_small
    memory "${params.m_mem} GB"
    cpus params.mtp_cores
    
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





