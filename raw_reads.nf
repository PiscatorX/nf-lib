#!/usr/bin/env nextflow




process  ubam2fastq{

    publishDir "$launchDir/${params.WD}/ubam2fastq/"
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

    publishDir "$launchDir/${params.WD}/fastq_SE/"
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

        "$launchDir/${params.WD}/fastq_SE/" 
	input:
		path fastqc_data

	output:
		path "multiqc*"


	"""
	
           multiqc .


	"""
	 
}
