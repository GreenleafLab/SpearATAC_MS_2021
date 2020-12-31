# High-throughput single-cell chromatin accessibility CRISPR screens enable unbiased identification of regulatory networks in cancer (Pierce SE*, Granja JM*, et al. 2020)

## **Link** : https://www.biorxiv.org/content/10.1101/2020.11.02.364265v1

## Please cite : Pierce SE et al., High-throughput single-cell chromatin accessibility CRISPR screens enable unbiased identification of regulatory networks in cancer. bioRxiv (2020) <br/>

![](Images/Figure1.png)

# Brief Descriptions of Analysis Scripts

## SPEAR-ATAC Analyses (See Scripts Folder)

**Download-Test-SpearATAC-Data.R - Code for downloading SpearATAC test data.** 

**SpearATAC-Align-sgRNA.R** - Code for analysis of SpearATAC Screen. 

**Analysis-of-Spear-ATAC-Screen.R** - Code for analysis of SpearATAC Screen. 

**SpearATAC-Functions.R** - Helper Functions for SpearATAC Analysis. 

# Test Data is From K562 Large Screen

## Data Files

**data/K562_R1.fragments.tsv.gz** - 10x Genomics Fragment Files <br/>
**data/K562_R2.fragments.tsv.gz** - 10x Genomics Fragment Files <br/>

**data/K562_R1.singlecell.csv** - 10x Genomics Summary Stats containing cells passing filters <br/>
**data/K562_R1.singlecell.csv** - 10x Genomics Summary Stats containing cells passing filters <br/>

**data/K562_R1.sgRNA.rds** - SpearATAC sgRNA Files containing Barcode Alignments <br/>
**data/K562_R2.sgRNA.rds** - SpearATAC sgRNA Files containing Barcode Alignments <br/>

**Vierstra-Human-Motifs.rds** - Motifs from Vierstra et al Nature 2020