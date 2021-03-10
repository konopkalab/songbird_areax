# Single-nucleus transcriptomics in Area X of male zebra finches (_Taeniopygia guttata_)
Resources for single nucleus RNA-sequencing data analysis presented in:

Xiao, L., Merullo, D. P., Cao, M., Co, M., Kulkarni, A., Konopka, G., & Roberts, T. F. (2020). Expression of FoxP2 in the basal ganglia regulates vocal motor sequences in the adult songbird. _BioRxiv_. https://doi.org/10.1101/2020.03.14.991042

# Directories

- MCX20: Raw data for the scrambled control group (CS-shScr) to be loaded with `Read10X` in Seurat. 

- MCX21: Raw data for the experimental FoxP2 knockdown group (CS-shFoxP2) to be loaded with `Read10X` in Seurat. 

- data:<br>
  - 01-dataclustering.md: Markdown file of R code to create Seurat objects from raw data. Processed Seurat objects can be downloaded [at this link](https://cloud.biohpc.swmed.edu/index.php/s/nLicEtkmjGGmRF8).
  
  - 02-code: R code to produce all figures (main and supplementary) in the accompanying study. 
  
  - All other files: Processed output from analyses; to be loaded when needed when running code. 
