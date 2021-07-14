#!/usr/bin/env nextflow




process  ubam2fastq{

    publishDir "$launchDir/${params.WD}/ubam2fastq/", mode: 'copy'
    echo true
    input:
	path uBAM 

   output:
       file fastq 

    script:
      fastq =   uBAM.getName().replace('bam','fastq') 

    """

       bedtools \
       		bamtofastq \
       		-i ${uBAM} \
		-fq ${fastq}

    """

}



process fastqc_SE{

    publishDir "$launchDir/${params.WD}/fastq_SE/", mode: 'move'
    input:
	path SE_reads

    output:
       file "${SE_reads.baseName}*"

    """
      
      fastqc \
      	     --extract \
    	     -f fastq \
    	     -t 2 \
             ${SE_reads}

    """
    
}



process multiqc{

        publishDir "$launchDir/${params.WD}/fastq_SE/", mode: 'move' 
	input:
		path fastqc_data

	output:
		path "multiqc*"


	"""
	
           multiqc .


	"""
	 
}




process trimmomatic_SE{

    memory "${params.m_mem} GB"
    cpus  params.mtp_cores
    publishDir path: "$launchDir/${params.WD}/trimmomatic", mode: 'copy'
    
    input:
	path SE_read 

    output:
        path ("trim_${SE_read}"), emit: TrimmedRead
        path ("${read_basename}.*")  

    script:
	read_basename =  SE_read.getSimpleName()
    
"""
	
   trimmomatic SE \
	${SE_read} \
        trim_${SE_read} \
	-threads ${params.mtp_cores} \
	-phred33\
	-trimlog ${read_basename}.log \
        -summary ${read_basename}.summary \
	LEADING:10 \
	TRAILING:10 \
	SLIDINGWINDOW:25:10 \
	MINLEN:50  
   
"""	

}
