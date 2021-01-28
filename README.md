# Single-nucleus transcriptomics in Area X of male zebra finches (_Taeniopygia guttata_)
Resources for single nuclei RNA-seq data analysis presented in:

Xiao, L., Merullo, D. P., Cao, M., Co, M., Kulkarni, A., Konopka, G., & Roberts, T. F. (2020). Expression of FoxP2 in the basal ganglia regulates vocal motor sequences in the adult songbird. _BioRxiv_. https://doi.org/10.1101/2020.03.14.991042

# Directories

- MCX20: Raw data for the scrambled control group (CS-shScr) to be loaded with `Read10X` in Seurat. 

- MCX21: Raw data for the experimental FoxP2 knockdown group (CS-shFoxP2) to be loaded with `Read10X` in Seurat. 

- data:<br>
  - dataclustering.md: Markdown file of R code for analysis.

  - MCX20_C_counttable.txt: Table of information exported from Seurat object used for analysis of the CS-shScr dataset.
    - Cell: Individual cell barcode<br>
    - UMIs: Number of unique molecular identifiers<br>
    - Genes: Number of genes<br>
    - Percent Mito Genes: Percentage of genes that are mitochondrial genes<br>
    - Cluster: Assigned cluster from Seurat analysis<br>
    - All other columns: mRNA transcript counts for each gene
  
  - MCX_COMBINED_C_seurat_FOXP2_markers.csv: Table of differentially expressed genes in FoxP2+ cells between CS-shScr and CS-shFoxP2 groups. [See the Seurat vignette on differential expression testing for column meanings](https://satijalab.org/seurat/v3.0/de_vignette.html):  
    > - p_val: p_val (unadjusted)<br>
    > - avg_logFC: log fold-change of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group<br>
    > - pct.1: The percentage of cells where the feature is detected in the first group<br>
    > - pct.2: The percentage of cells where the feature is detected in the second group<br>
    > - p_val_adj : Adjusted p-value, based on bonferroni correction using all features in the dataset.
