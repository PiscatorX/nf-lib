#!/usr/bin/env nextflow

process  ubam2fastq{

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
      
       head ${fastq}

    """


}
