process catx{

    memory "${params.m_mem} GB"
    input:
        path some_files 

    output:
	path "combined.files", emit: combined  		


"""

   cat  *  >  combined.files  

"""
}



process merge_tsv{

    memory "${params.l_mem} GB"
    publishDir path: "${params.WD}/${outdir}", mode: 'copy'
    input:
	path infoseq
	val outdir

    output:
	path "merged_df.tsv"

"""  
 #!/usr/bin/env python
import pandas as pd
import glob

globfiles = glob.glob("*")

merged_df = pd.concat([ pd.read_csv(tsv, sep ="\t") for tsv  in  globfiles ])

print(merged_df.describe())


merged_df.to_csv("merged_df.tsv")
			 	     
"""

}



