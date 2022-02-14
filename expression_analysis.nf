

process transcript_est_alignfree{

    cpus params.ltp_cores 
    memory "${params.l_mem} GB"    
    publishDir path: "$params.WD/trinity-${est_method}-transcript-est",  mode: 'copy'
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
       --thread_count ${params.mtp_cores} 

"""

}




process abundance_estimates_to_matrix{

//https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats

    cpus params.ltp_cores 
    memory "${params.l_mem} GB"    
    publishDir path: "$params.WD/trinity-${est_method}-abund-matrix",  mode: 'move'    
    input:
        path quant_dir
	path denovo_ref
	val quant_name
	val est_method

    output:
        path "${est_method}*", emit: matrix_data
	path "ExN50.stats", emit: ExN50

"""

 find -L . -iname "${quant_name}" -print  >  quant_files

 abundance_estimates_to_matrix.pl \
    --est_method salmon \
    --quant_files quant_files \
    --gene_trans_map none \
    --name_sample_by_basedir      

 contig_ExN50_statistic.pl \
     ${est_method}.isoform.TMM.EXPR.matrix ${denovo_ref} > ExN50.stats

"""
//cat transcripts.TMM.EXPR.matrix.E-inputs |  egrep -v ^\# | awk '$1 <= 90' | wc -l
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
//https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification


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
        --gcBias \
	--no-version-check \
        --alignments ${bam} \
        --targets ${fasta_reference} \
	--writeUnmappedNames \
	--output ${bam.baseName} \
	--threads ${params.htp_cores}

    mv ${bam.baseName}/quant.sf  ${bam.baseName}.sf    

"""	
}
