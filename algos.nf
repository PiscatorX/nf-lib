#!/usr/bin/env nextflow



process cd_hit_est {

     cpus params.htp_cores
     memory "${params.m_mem} GB"
     publishDir path: "${params.WD}/CD-Hit/${tag}/${contig_basename}", mode: 'copy'
     input:
          val perc
          path fasta_seqs
          val tag
       
    output:
         path "${contig_basename}.cd_hits", emit: cd_hits 
         path "${contig_basename}.cd_hits.clstr", emit: cdhit_clusters
         


    script:
 	contig_basename = "${fasta_seqs.baseName}"

"""
    cd-hit-est \
        -i ${fasta_seqs} \
	-c ${perc} \
	-T ${params.htp_cores} \
	-M 0 \
	-d 0 \
	-r 0 \
	-p 1 \
	-g 1 \
	-o ${contig_basename}.cd_hits 

"""

}


process cd_hit_consensus{

     cpus params.ltp_cores
     memory "${params.m_mem} GB"
     publishDir path: "${params.WD}/CD-Hit/${tag}/", mode: 'copy'
     input:
          path cd_hits
          path cd_hits_clstr
	  val min_size
          val  tag
       
    output:
          path "CDhit_consensus.fasta", emit: cd_hit_consensus 

    script:
 	cdhits_basename = "${cd_hits.baseName}"
	
"""
    make_multi_seq.pl \
        ${cd_hits} \
        ${cd_hits_clstr} \
        SeqDir \
        ${min_size}


   for file in \$(ls SeqDir/*)
     do 
       cat  \${file} >>  CDhit_consensus.fasta 
     done


"""
//To run this script correctly,
//"-d 0" option should be used in the cd-hit run and it is better to use "-g 1"

}



