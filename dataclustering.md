# DATA CLUSTERING FOR SCRAMBLED GROUP

This script clusters the data from the Scrambled group. The Seurat pipeline (v3.0) is used to filter, normalize, scale, cluster, and identify gene markers. 

**LOAD DATA AND CREATE SEURAT OBJECT**

Unfiltered data are imported after processing through the Cellranger pipeline.

```
# LOAD DATA
# INTERNAL NAME FOR SCRAMBLED GROUP IS "MCX20C"
MCX20C <- 
  Read10X("MCX20C")

# CREATE SEURAT OBJECT
MCX20C_seurat <-
  CreateSeuratObject(counts = MCX20C)
```

**FILTERING THE DATA**

Filter data based on thresholds for the number of UMIs and the percent mitochondrial genes. 

```
# ADD MITOCHONDRIAL DATA
MCX20C_seurat_mito_genes <- 
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

MCX20C_seurat_percent_mito <-
  Matrix::colSums(
    GetAssayData(MCX20C_seurat)[MCX20C_seurat_mito_genes, ])/Matrix::colSums(
      GetAssayData(MCX20C_seurat)
      )

MCX20C_seurat$percent.mito <-
  MCX20C_seurat_percent_mito

# EXAMINE DIAGNOSTIC VIOLIN PLOTS
VlnPlot(
  MCX20C_seurat,
  features = c(
  "nFeature_RNA", 
  "nCount_RNA", 
  "percent.mito"
  )
)

# EXAMINE DIAGNOSTIC GENE PLOTS
MCX20C_plot1 <-
  FeatureScatter(
    MCX20C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "percent.mito"
  )

MCX20C_plot2 <- 
  FeatureScatter(
    MCX20C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

CombinePlots(plots = list(MCX20C_plot1, MCX20C_plot2))

# FILTER CELLS BY NCOUNT AND PERCENT MITOCHONDRIAL
MCX20C_seurat <-
  subset(
    MCX20C_seurat,
    subset = 
    	nCount_RNA > -Inf & 
    	nCount_RNA < 5000 &
    	percent.mito < 0.05
  )
```

**NORMALIZE AND SCALE THE DATA**

Normalize and scale the data using the Seurat pipeline. 

```
# NORMALIZE DATA
MCX20C_seurat <-
  NormalizeData(
    MCX20C_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

# IDENTIFY THE TOP 2000 VARIABLE GENES
MCX20C_seurat <-
  FindVariableFeatures(MCX20C_seurat)

# PLOT VARIABLE GENES
VariableFeaturePlot(MCX20C_seurat)

# SCALE DATA, REGRESSING ON THE NUMBER OF UMIS AND PCT MITOCHONDRIAL GENES
MCX20C_seurat <-
  ScaleData(
    MCX20C_seurat,
    vars.to.regress = c("nCount_RNA", "percent.mito"),
    model.use = "linear"
  )
```

**PCA ANALYSIS**

The top contributing PCs are identified to be used in downstream analyses.

```
# RUN PCA WITH TOP 2000 VARIABLE GENES
MCX20C_seurat <-
  RunPCA(
    MCX20C_seurat,
    features = VariableFeatures(object = MCX20C_seurat),
    ndims.pint = 1:5
  )

# VISUALIZE CONTRIBUTIONS OF PCS WITH ELBOW PLOT
ElbowPlot(MCX20C_seurat)
```

**DATA CLUSTERING**

A cluster plot is produced following the Seurat pipeline. Clusters are assigned with the Louvain algorithm and then cells are plotted in UMAP space.

```
# CONSTRUCT A SHARED NEAREST NEIGHBOR (SNN) GRAPH WITH IDENTIFIED PCS
MCX20C_seurat <- 
  FindNeighbors(
    MCX20C_seurat,
    dims = 1:18
  )

# CLUSTER CELLS WITH THE LOUVAIN ALGORITHM 
MCX20C_seurat <-
  FindClusters(
    MCX20C_seurat,
    algorithm = 1,
    resolution = 1.2,
  )

# RUN UMAP
MCX20C_seurat <-
  RunUMAP(
    MCX20C_seurat,
    dims = 1:18
  )

# CREATE PLOT OF UMAP OUTPUT
MCX20C_seurat_clusterplot <-
  UMAPPlot(
    MCX20C_seurat,
    label = TRUE
    )

# SAVE UMAP PLOT
save(
  MCX20C_seurat,
  file = "MCX20C_seurat.RData"
)
```

**CLUSTER MARKERS**

Identify differentially expressed gene markers for each cluster.

```
#RUN DEG TEST
MCX20C_markers <- 
  FindAllMarkers(
    MCX20C_seurat, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )

#SAVE AS DATA.TABLE
write.table(
  MCX20C_markers,
  file = "MCX20C_markers.txt"
  ),
  row.names = FALSE
)
```



# DATA CLUSTERING FOR INTEGRATED DATA

First, data are loaded from the Scrambled and Knockdown groups following the process described above. 

**Scrambled Group**

```
# LOAD DATA
# INTERNAL NAME FOR SCRAMBLED GROUP IS "MCX20C"
MCX20C <- 
  Read10X("MCX20C")
  
MCX20C_seurat <-
  CreateSeuratObject(
    counts = MCX20C, 
    project = "MCX20C_seurat",
  )

MCX20C_seurat$dataset <- 
  "MCX20C"

MCX20C_seurat_mito_genes <- 
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

MCX20C_seurat_percent_mito <-
  Matrix::colSums(
    GetAssayData(MCX20C_seurat)[MCX20C_seurat_mito_genes, ])/Matrix::colSums(
      GetAssayData(MCX20C_seurat))

MCX20C_seurat$percent.mito <-
  MCX20C_seurat_percent_mito

MCX20C_seurat <- 
  subset(
    MCX20C_seurat, 
    subset = nCount_RNA < 5000 & percent.mito < 0.05
  )

MCX20C_seurat <-
  NormalizeData(
    MCX20C_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

MCX20C_seurat <-
  FindVariableFeatures(MCX20C_seurat)
```

**Knockdown Group**

```
# LOAD DATA
# INTERNAL NAME FOR KNOCKDOWN GROUP IS "MCX21C"
MCX21C <- 
  Read10X("MCX21C")

MCX21C_seurat <-
  CreateSeuratObject(
    counts = MCX21C, 
    project = "MCX21C_seurat",
  )

MCX21C_seurat$dataset <- 
  "MCX21C"

MCX21C_seurat_mito_genes <- 
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

MCX21C_seurat_percent_mito <-
  Matrix::colSums(
    GetAssayData(MCX21C_seurat)[MCX21C_seurat_mito_genes, ])/Matrix::colSums(
      GetAssayData(MCX21C_seurat))

MCX21C_seurat$percent.mito <-
  MCX21C_seurat_percent_mito

MCX21C_seurat <- 
  subset(
    MCX21C_seurat, 
    subset = nCount_RNA < 5000 & percent.mito < 0.05
  )

MCX21C_seurat <-
  NormalizeData(
    MCX21C_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

MCX21C_seurat <-
  FindVariableFeatures(MCX21C_seurat)
```

**INTEGRATION**

The Scrambled and Knockdown datasets are integrated using the Seurat pipeline. The Integrated data assay is then set for all downstream analyses. 

```
# INTERNAL NAME FOR INTEGRATED DATA IS "MCXC"
# FIND INTEGRATION ANCHORS
MCXC_seurat_anchors <-
  FindIntegrationAnchors(
    object.list = list(
      MCX20C_seurat,
      MCX21C_seurat
    ),
    dims = 1:20
  )

# CREATE LIST OF ALL GENES TO BE INTEGRATED
to_integrate <- 
  Reduce(
    intersect,
    lapply(
      MCXC_seurat.anchors@object.list,
      rownames
    )
  )

# INTEGRATE DATA USING ALL GENES
MCXC_seurat <-
  IntegrateData(
    anchorset = MCXC_seurat_anchors,
    features.to.integrate = to_integrate,
    dims = 1:20
  )

# SET THE DEFAULT ASSAY TO BE INTEGRATED FOR DOWNSTREAM STEPS
DefaultAssay(MCXC_seurat) <-
  "integrated"
```

**SCALING, PCA, CLUSTERING, AND UMAP**

These steps are performed as described in the beginning section on clustering for the Scrambled group.
```
# SCALE DATA
MCXC_seurat <-
  ScaleData(
    MCXC_seurat,
    vars.to.regress = c(
      "nCount_RNA", "percent.mito"
    ),
    model.use ="linear"
  )

# RUN PCA
MCXC_seurat <-
  RunPCA(
    MCXC_seurat,
    features = VariableFeatures(
      object = MCXC_seurat
    ),
    ndims.pint = 1:5
  )

# PC ELBOW PLOT
ElbowPlot(MCXC_seurat)

# CLUSTER PLOT
MCXC_seurat <- 
  FindNeighbors(
    MCXC_seurat,
    dims = 1:16
  )

MCXC_seurat <-
  FindClusters(
    MCXC_seurat,
    resolution = 0.8,
  )

MCXC_seurat <-
  RunUMAP(
    MCXC_seurat,
    dims = 1:16
  )

MCXC_seurat_clusterplot <-
  UMAPPlot(
    MCXC_seurat,
    label = TRUE)

save(
  MCXC_seurat,
  file = "MCXC_seurat.RData"
)
```
