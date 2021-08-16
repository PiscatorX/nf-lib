process Trinity_SE{

    cpus params.htp_cores
    memory "${params.h_mem} GB"
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

