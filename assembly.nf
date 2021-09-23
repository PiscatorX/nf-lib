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
     

 		   
"""        

     Trinity\
	 --seqType fq \
	 --samples_file ${params.sample_file} \
	 --max_memory ${params.h_mem}G \
	 --CPU ${params.htp_cores} \
	 --output Trinity \
	 --verbose
     
    TrinityStats.pl Trinity/Trinity.fasta > Trinity/Assembly.stats
   
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

    publishDir path: "$params.WD",  mode: 'move'
    params.scratch_small
    cpus params.htp_cores
    input:
        path SE_reads
	path fasta_reference
	path genome_ref

    output:
        path "transrate_results"
	stdout emit: transrate_stdout
	

script:
Left_reads = SE_reads.collect{it }.join(', ')

"""

   transrate \
       --left ${Left_reads} \
       --assembly ${fasta_reference} \
       --threads ${params.htp_cores} \
       --reference ${genome_ref}

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
      -c TRANSCRIPTS [TRANSCRIPTS ...]] [-psl ALIGNMENT [ALIGNMENT ...]]
                        [-sam READS_ALIGNMENT] [-1 LEFT_READS] [-2 RIGHT_READS] [-s SINGLE_READS] [--gmap_index GMAP_INDEX] [-o OUTPUT_DIR] [--test] [-d] [-t THREADS]
                        [-l LABELS [LABELS ...]] [-ss] [--min_alignment MIN_ALIGNMENT] [--no_plots] [--blat] [--gene_mark] [--meta] [--lower_threshold LOWER_THRESHOLD]
                        [--upper_threshold UPPER_THRESHOLD] [--disable_infer_genes] [--disable_infer_transcripts] [--busco BUSCO] [--prokaryote]
 
       
       
"""

}

 // /opt/rnaQUAST.py [-h] [-r REFERENCE [REFERENCE ...]] [--gtf GTF [GTF ...]] [--gene_db GENE_DB] [-c TRANSCRIPTS [TRANSCRIPTS ...]] [-psl ALIGNMENT [ALIGNMENT ...]]
 //                        [-sam READS_ALIGNMENT] [-1 LEFT_READS] [-2 RIGHT_READS] [-s SINGLE_READS] [--gmap_index GMAP_INDEX] [-o OUTPUT_DIR] [--test] [-d] [-t THREADS]
 //                        [-l LABELS [LABELS ...]] [-ss] [--min_alignment MIN_ALIGNMENT] [--no_plots] [--blat] [--gene_mark] [--meta] [--lower_threshold LOWER_THRESHOLD]
 //                        [--upper_threshold UPPER_THRESHOLD] [--disable_infer_genes] [--disable_infer_transcripts] [--busco BUSCO] [--prokaryote]




