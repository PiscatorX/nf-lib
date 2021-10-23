process Trinity_SE{

    cpus params.htp_cores
    memory "${params.h_mem} GB"
    scratch params.scratch_large
    publishDir "$params.WD",  mode: 'move'
    
    input:
	path trimmed_fastq_SE
	path sample_file
	
    output:
         path "Trinity"
	 path "Trinity/Trinity.fasta", emit: assembly    

 		   
"""        

     Trinity\
	 --seqType fq \
	 --samples_file ${params.sample_file} \
	 --max_memory ${params.h_mem}G \
	 --CPU ${params.htp_cores} \
	 --output Trinity \
	 --verbose

"""

}



process assembly_stats{

    cpus params.ltp_cores
    memory "${params.l_mem} GB"
    publishDir "$params.WD/Trinity",  mode: 'move'
    
    input:
	path assembly
	
    output:
         path "Assembly.stats"
	 
"""

   TrinityStats.pl ${assembly} > Assembly.stats

"""

}




process analyze_blastTophits{

    cpus params.ltp_cores
    scratch params.scratch_small
    memory "${params.l_mem} GB"    
    publishDir path: "$params.WD/Assembly",  mode: 'move'
    input:
        path blastout_fmt6
	path denovo_ref
	path prot_reference
	
    output:
        path "${blastout_fmt6}.*" 

	    
"""

    analyze_blastPlus_topHit_coverage.pl \
        ${blastout_fmt6} \
        ${denovo_ref} \
        ${prot_reference}  >  ${blastout_fmt6}.coverage


    $TRINITY_HOME/util/misc/blast_outfmt6_group_segments.pl \
        ${blastout_fmt6} \
        ${denovo_ref} \
        ${prot_reference}  >  ${blastout_fmt6}.grouped

    $TRINITY_HOME/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl \
         ${blastout_fmt6}.grouped > ${blastout_fmt6}.grouped.coverage

	      
"""

}



process transrate{

    echo true
    publishDir path: "$params.WD/Transrate",  mode: 'move'
    params.scratch_small
    cpus params.htp_cores
    input:
        path SE_reads
	val fasta_reference
	val genome_ref

    output:
        path "transrate_results"
	path "transrate_stdout"
	

script:
Left_reads = SE_reads.collect{it }.join(', ')

"""

    transrate \
	--left ${Left_reads} \
	--assembly ${fasta_reference} \
	--threads ${params.htp_cores} \
	--reference ${genome_ref}  > transrate_stdout

"""

}




process rsem_eval_est{

    publishDir path: "$params.WD/detonate_results",  mode: 'move'
    cpus params.htp_cores
    input:
        path SE_reads
	path fasta_reference
	path genome_ref

	
script:
Left_reads = SE_reads.collect{it }.join(', ')
(readlen_mean, readlen_stddev) = readlen_stats.split("\t")
"""

    rsem-eval-estimate-transcript-length-distribution \
        ${fasta_reference} \
        ${readlen_men} \
        parameter_file   

"""

}




process rsem_eval{

    publishDir path: "$params.WD/detonate_results",  mode: 'move'
    cpus params.htp_cores

    input:
        path SE_reads
	path fasta_reference
	path genome_ref
	val  readlen_stats
	
script:
SE_reads = SE_reads.collect{it }.join(', ')
(readlen_mean, readlen_stddev) = readlen_stats.split("\t")
"""

   rsem-eval-estimate-transcript-length-distribution \
       ${fasta_reference} \
       ${fasta_reference.baseName}.txt

   rsem-eval-calculate-score \
       --num-threads  ${params.htp_cores} \
       ${SE_reads} \
       ${fasta_reference} \
       ${fasta_reference.baseName} \
       ${readlen_mean} 
       
       
"""

}




process rnaQuast{

    echo true
    publishDir path: "$params.WD/rnaQuast",  mode: 'move'
    cpus params.htp_cores
    input:
        path SE_reads
	path fasta_reference
	path genome_ref
	val  readlen_stats
	
script:
SE_reads = SE_reads.collect{it }.join(', ')
(readlen_mean, readlen_stddev) = readlen_stats.split("\t")
"""

    rnaQUAST.py \
        --transcripts TRANSCRIPTS \
        --reference REFERENCE \
        --gtf GENE_COORDINATES
   
       
"""

}

 


