



process  ubam2fastq{

    publishDir "${params.WD}/ubam2fastq/", mode: 'copy'
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

    publishDir "${params.WD}/fastq_SE/${outdir}"
    cpus params.ltp_cores
    input:
    	path SE_reads
	val outdir
	
	
    output:
        path "${fastqc_basename}*"


    script:
        fastqc_basename = SE_reads.getSimpleName()

    """
        
         fastqc \
    	     -f fastq \
    	     --threads ${params.ltp_cores} \
             ${SE_reads}

    """
    
}



process multiqc{

        publishDir "${params.WD}/multiqc/", mode: 'move' 
	input:
	     path fastqc_html
	     val outdir

	output:
             path "${outdir}*"


	"""
            multiqc  . \
            --filename ${outdir} 

	"""
}



process fix_ReadName{

        echo true
	input:
	     path SE_read

	output:
	     path SE_read
shell:
"""
     
     sed -i -e '1~4s/\$/\\/1/g' ${SE_read} 
     

"""

}



process trimmomatic_SE{

    memory "${params.m_mem} GB"
    cpus  params.mtp_cores
    publishDir path: "${params.WD}/trimmomatic", mode: 'copy'
    
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


