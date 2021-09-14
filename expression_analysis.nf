

process transcript_est_alignfree{

    cpus params.ltp_cores 
    memory "${params.l_mem} GB"    
    publishDir path: "$params.WD/${est_method}",  mode: 'move'
    input:
        path trimmed_fastq_SE
        path denovo_ref
	path sample_file
	val readlen_stats
	val est_method

    output:
        path "*rep*", emit: quant_dir
  	path "*idx"   
	
script:
(readlen_mean, readlen_stddev) = readlen_stats.split("\t")
"""
    
   align_and_estimate_abundance.pl \
       --prep_reference \
       --est_method ${est_method} \
       --transcripts ${denovo_ref} \
       --thread_count ${params.mtp_cores} 

   align_and_estimate_abundance.pl \
       --samples_file ${sample_file} \
       --output_dir ${est_method} \
       --est_method ${est_method} \
       --transcripts ${denovo_ref} \
       --seqType fa \
       --trinity_mode \
       --thread_count ${params.mtp_cores} 

"""

}


   






process salmon_index{
    
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    storeDir "${DB_REF}/Salmon/"

    input:
	path reference

    output:
        path "${salmon_index}*", emit: salmon_index

script:
salmon_index = reference.getSimpleName()
"""

    salmon \
        index \
        --no-version-check \
        -t  ${reference}  \
	-i  ${salmon_index} \
	--type quasi \
	-p  ${params.htp_cores} 

"""


}



process salmon_quant{
    
    cpus params.mtp_cores
    memory "${params.h_mem} GB"
    publishDir path: "$params.WD/Salmon/", pattern: "${bam.baseName}", mode: 'copy'
    publishDir path: "$params.WD/Salmon-Quant/", pattern: "${bam.baseName}.sf", mode: 'copy'
 
    input:
        path bam
	path salmon_index
	path fasta_reference
	
    output:
        path "${bam.baseName}"
	path "${bam.baseName}.sf", emit: salmon_quant

"""

    salmon \
        quant \
        --libType A \
	--no-version-check \
        --alignments ${bam} \
        --targets ${fasta_reference} \
	--writeUnmappedNames \
	--output ${bam.baseName} \
	--threads ${params.htp_cores}

    mv ${bam.baseName}/quant.sf  ${bam.baseName}.sf    

"""

	
}

