**DATA CLUSTERING FOR SCRAMBLED GROUP**

This script clusters the data from the Scrambled group. The Seurat pipeline (v3.0) is used to filter, normalize, scale, cluster, and identify gene markers. 

**LOAD DATA AND CREATE SEURAT OBJECT**

Unfiltered data are imported after processing through the Cellranger pipeline.

```
#LOAD DATA
#INTERNAL NAME FOR SCRAMBLED GROUP IS "MCX20_C"
MCX20_C <- 
  Read10X("MCX20C")

#CREATE SEURAT OBJECT
MCX20C_seurat <-
  CreateSeuratObject(counts = MCX20C)
```

**FILTERING THE DATA**

Filter data based on thresholds for the number of UMIs and the percent mitochondrial genes. 

```
#ADD MITOCHONDRIAL DATA
MCX20_C_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

MCX20_C_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    MCX20_C_seurat,
    features = MCX20_C_seurat_mito_genes
  )

#DIAGNOSTIC VIOLIN PLOTS
VlnPlot(
  MCX20_C_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)

#DIAGNOSTIC GENE PLOTS
MCX20_C_seurat_plot1 <-
  FeatureScatter(
    MCX20_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "percent.mito"
  )

MCX20_C_seurat_plot2 <- 
  FeatureScatter(
    MCX20_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

CombinePlots(
  plots = list(
    MCX20_C_seurat_plot1,
    MCX20_C_seurat_plot2
    )
  )

#FILTER CELLS BY NCOUNT AND PERCENT MITOCHONDRIAL
MCX20_C_seurat <-
  subset(
    MCX20_C_seurat,
    subset = nCount_RNA > -Inf & nCount_RNA < 10000 & percent.mito < 5
  )
```

**NORMALIZE AND SCALE THE DATA**

Normalize and scale the data using the Seurat pipeline. 

```
#NORMALIZE DATA
MCX20_C_seurat <-
  NormalizeData(
    MCX20_C_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#IDENTIFY THE TOP 2000 VARIABLE GENES
MCX20_C_seurat <-
  FindVariableFeatures(MCX20_C_seurat)

#PLOT VARIABLE GENES
VariableFeaturePlot(MCX20_C_seurat)

#SCALE DATA
MCX20_C_seurat <-
  ScaleData(
    MCX20_C_seurat,
    vars.to.regress = c("nCount_RNA", "percent.mito"),
    model.use = "linear"
  )
```

**PCA ANALYSIS**

The top contributing PCs are identified to be used in downstream analyses.

```
#RUN PCA
MCX20_C_seurat <-
  RunPCA(
    MCX20_C_seurat,
    features = VariableFeatures(object = MCX20_C_seurat),
    ndims.pint = 1:5
  )

#PC ELBOW PLOT
ElbowPlot(MCX20_C_seurat)
```

**DATA CLUSTERING**

A cluster plot is produced following the Seurat pipeline. Clusters are assigned with the Louvain algorithm and then cells are plotted in UMAP space.

```
#CONSTRUCT A SHARED NEAREST NEIGHBOR (SNN) GRAPH WITH IDENTIFIED PCS
MCX20_C_seurat <- 
  FindNeighbors(
    MCX20_C_seurat,
    dims = 1:15
  )

#CLUSTER CELLS WITH THE LOUVAIN ALGORITHM 
MCX20_C_seurat <-
  FindClusters(
    MCX20_C_seurat,
    resolution = 0.6,
  )
  
#RUN UMAP
MCX20_C_seurat <-
  RunUMAP(
    MCX20_C_seurat,
    dims = 1:15
  )

#CREATE PLOT OF UMAP OUTPUT
  UMAPPlot(
    MCX20_C_seurat,
    label = TRUE
  )

#SAVE UMAP PLOT
save(
  MCX20_C_seurat,
  file = "MCX20_C_seurat.RData"
)
```

**CLUSTER MARKERS**

Identify differentially expressed gene markers for each cluster.

```
#RUN DEG TEST
MCX20_C_seurat_markers <- 
  FindAllMarkers(MCX20_C_seurat)

#SAVE AS DATA.TABLE
write.table(
  MCX20_C_seurat_markers,
  file = "MCX20_C_seurat_markers.txt"
  ),
  row.names = FALSE
)
```



**DATA CLUSTERING FOR INTEGRATED DATA**

First, data are loaded from the Scrambled and Knockdown groups following the process described above. 

**Scrambled Group**

```
#LOAD DATA
#INTERNAL NAME FOR SCRAMBLED GROUP IS "MCX20C"
MCX20_C <- 
  Read10X("MCX20C")
  
#SEURAT OBJECT
MCX20_C_seurat <-
  CreateSeuratObject(
    counts = MCX20_C,
    project = "MCX20"
  )

#ADD MITOCHONDRIAL DATA
MCX20_C_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

MCX20_C_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    MCX20_C_seurat,
    features = MCX20_C_seurat_mito_genes
  )

#DIAGNOSTIC VIOLIN PLOTS
VlnPlot(
  MCX20_C_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)

#DIAGNOSTIC GENE PLOTS
plot1 <-
  FeatureScatter(
    MCX20_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "percent.mito"
  )

plot2 <- 
  FeatureScatter(
    MCX20_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

CombinePlots(plots = list(plot1, plot2))

#FILTER CELLS BY NCOUNT AND PERCENT MITOCHONDRIAL
MCX20_C_seurat <-
  subset(
    MCX20_C_seurat,
    subset = nCount_RNA > -Inf & nCount_RNA < 10000 & percent.mito < 5
  )

#NORMALIZE DATA
MCX20_C_seurat <-
  NormalizeData(
    MCX20_C_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
MCX20_C_seurat <-
  FindVariableFeatures(MCX20_C_seurat)
```

**Knockdown Group**

```
#LOAD DATA
#INTERNAL NAME FOR KNOCKDOWN GROUP IS "MCX21_C"
MCX21_C <- 
  Read10X("MCX21C")

#SEURAT OBJECT
MCX21_C_seurat <-
  CreateSeuratObject(
    counts = MCX21_C,
    project = "MCX21"
  )

#ADD MITOCHONDRIAL DATA
MCX21_C_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

MCX21_C_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    MCX21_C_seurat,
    features = MCX21_C_seurat_mito_genes
  )

#DIAGNOSTIC VIOLIN PLOTS
VlnPlot(
  MCX21_C_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)

#DIAGNOSTIC GENE PLOTS
plot1 <-
  FeatureScatter(
    MCX21_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "percent.mito"
  )

plot2 <- 
  FeatureScatter(
    MCX21_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

CombinePlots(plots = list(plot1, plot2))

#FILTER CELLS BY NCOUNT AND PERCENT MITOCHONDRIAL
MCX21_C_seurat <-
  subset(
    MCX21_C_seurat,
    subset = nCount_RNA > -Inf & nCount_RNA < 10000 & percent.mito < 5
  )

#NORMALIZE DATA
MCX21_C_seurat <-
  NormalizeData(
    MCX21_C_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
MCX21_C_seurat <-
  FindVariableFeatures(MCX21_C_seurat)
```

**INTEGRATION**

The Scrambled and Knockdown datasets are integrated using the Seurat pipeline. The Integrated data assay is then set for all downstream analyses. 

```
#INTERNAL NAME FOR INTEGRATED DATA IS "MCX_COMBINED_C"
#FIND INTEGRATION ANCHORS
MCX_COMBINED_C_seurat_anchors <-
  FindIntegrationAnchors(
    object.list = list(
      MCX20_C_seurat,
      MCX21_C_seurat
    ),
    dims = 1:20
  )

#INTEGRATE DATA
MCX_COMBINED_C_seurat <-
  IntegrateData(
    anchorset = MCX_COMBINED_C_seurat_anchors,
    dims = 1:20
  )

#SET THE DEFAULT ASSAY TO BE INTEGRATED FOR DOWNSTREAM STEPS
DefaultAssay(MCX_COMBINED_C_seurat) <-
  "integrated"
```

**SCALING, PCA, CLUSTERING, AND UMAP**

These steps are performed as described in the beginning section on clustering for the Scrambled group.
```
#SCALE DATA
MCX_COMBINED_C_seurat <-
  ScaleData(
    MCX_COMBINED_C_seurat,
    vars.to.regress = c(
      "nCount_RNA", "percent.mito"
    ),
    model.use ="linear"
  )

#RUN PCA
MCX_COMBINED_C_seurat <-
  RunPCA(
    MCX_COMBINED_C_seurat,
    features = VariableFeatures(
      object = MCX_COMBINED_C_seurat
    ),
    ndims.pint = 1:5
  )

#PC ELBOW PLOT
ElbowPlot(MCX_COMBINED_C_seurat)

#CLUSTER PLOT
MCX_COMBINED_C_seurat <- 
  FindNeighbors(
    MCX_COMBINED_C_seurat,
    dims = 1:15
  )

MCX_COMBINED_C_seurat <-
  FindClusters(
    MCX_COMBINED_C_seurat,
    resolution = 0.4,
  )

MCX_COMBINED_C_seurat <-
  RunUMAP(
    MCX_COMBINED_C_seurat,
    dims = 1:15
  )

MCX_COMBINED_C_seurat_clusterplot <-
  UMAPPlot(
    MCX_COMBINED_C_seurat,
    label = TRUE)

save(
  MCX_COMBINED_C_seurat,
  file = "MCX_COMBINED_C_seurat.RData"
)
```
