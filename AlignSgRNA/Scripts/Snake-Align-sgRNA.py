#Snakemake to run cell ranger count using slurm submit on sherlock!
import os
import sys
import glob
import yaml
import numpy as np
import pandas as pd

metadata_file = "metadata_align.txt"
if os.path.isfile(metadata_file):
    metadata = pd.read_table(metadata_file, index_col = False)
else:
    error = "metadata_align.txt does not exist! exiting..."
    sys.exit(error)

onsuccess:
    print("Snakemake finished successfully with no errors!")

sample_names = metadata.Name.tolist()
barcodes = "data/10x-Barcodes.txt"
sgRNA = "data/sgRNA-LargeScreen.txt"

rule all:
    input:
    	expand("output/Match_sgRNA/Chunks/{sample}-R1-match-sgRNA.rds", sample=sample_names),
    output:
        final = "Finished_Snakemake.txt"
    resources:
        mem_gb = 5,
        large_mem_gb = lambda wildcards, attempt: attempt * 5,    
    shell:
        "echo created by Jeffrey Granja. >> {output}"

rule Match_sgRNA_Read1:
    input:
        read1 = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample]["Read1"],
        index = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample]["Index"],
    output:
        "output/Match_sgRNA/Chunks/{sample}-R1-match-sgRNA.rds",
    threads: 2
    resources:
        mem_gb = 5,
    shell:
        "Rscript --vanilla Scripts/SpearATAC-Align-sgRNA.R "
        "--sgRNA {sgRNA} "
        "--barcodes {barcodes} "
        "--fastq_read {input.read1} "
        "--fastq_index {input.index} "
        "--output {output} "