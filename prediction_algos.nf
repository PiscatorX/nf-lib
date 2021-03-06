process transdecoder{

    publishDir path: "${params.WD}/transdcoder/", mode: 'copy'
    cpus  params.mtp_cores
    memory "${params.m_mem} GB"
    input:
	path denovo_ref
	path swissprot
	val swissprot_path
	val PfamA_hmm
	val evalue
	val min_len

    output:
        path "${denovo_ref}*.pep", emit: pep
	path "${denovo_ref}*"
        path "transdecoder*"
	path "blastp.outfmt6"
	path "pfam.domtblout"
	

"""

    TransDecoder.LongOrfs \
        -t ${denovo_ref} \
        -m ${min_len} \
        --output_dir transdecoder


    blastp \
        -query transdecoder/longest_orfs.pep  \
	-db ${swissprot_path}/${swissprot} \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue ${evalue} \
        -num_threads ${params.mtp_cores}  > blastp.outfmt6

    hmmscan \
        --cpu ${params.mtp_cores} \
        --domE ${evalue} \
        --domtblout pfam.domtblout \
        ${PfamA_hmm} \
        transdecoder/longest_orfs.pep

    TransDecoder.Predict \
        -t ${denovo_ref} \
        --output_dir transdecoder \
        --retain_pfam_hits pfam.domtblout \
        --retain_blastp_hits blastp.outfmt6
         
"""
}
