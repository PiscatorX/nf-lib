

process sam2bam{

    cpus params.mtp_cores 
    memory "${params.m_mem} GB"    
    publishDir path: "${params.WD}/samtools"
    input:
         path SAM


    output:
	path BAM

    script:
        BAM = SAM.getSimpleName() + '.sam'
	
"""

    samtools \
    	     view \
             -F 4 \
	     --threads ${params.mtp_cores} \
	     -b \
             -o ${BAM}

"""

}