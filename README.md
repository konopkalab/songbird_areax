# Single-nucleus transcriptomics in Area X of male zebra finches (_Taeniopygia guttata_)
Resources for single nuclei RNA-seq data analysis presented in:

Xiao, L., Merullo, D. P., Cao, M., Co, M., Kulkarni, A., Konopka, G., & Roberts, T. F. (2020). Expression of FoxP2 in the basal ganglia regulates vocal motor sequences in the adult songbird. _BioRxiv_. https://doi.org/10.1101/2020.03.14.991042

# Directories and Files

MCX20: Raw data for the scrambled control group (CS-shScr).

MCX21: Raw data for the experimental FoxP2 knockdown group (CS-shFoxP2).

dataclustering.md: Markdown file of R code for analysis.

MCX20_C_counttable.txt: Table of information exported from Seurat object used for analysis of the CS-shScr dataset.      
  --Cell: Individual cell barcode<br>
  --UMIs: Number of unique molecular identifiers<br>
  --Genes: Number of genes<br>
  --Percent Mito Genes: Percentage of genes that are mitochondrial genes<br>
  --Cluster: Assigned cluster from Seurat analysis
  --All other columns: mRNA transcript counts for each gene<br>
  
