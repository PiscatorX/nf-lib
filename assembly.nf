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