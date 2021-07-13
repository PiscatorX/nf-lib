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

    """

}



process fastqc_SE{

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

