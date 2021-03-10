library(tidyverse)
library(Seurat)
library(cowplot)
library(egg)
library(colorspace)
library(dendextend)
library(pvclust)
library(ggrepel)

#----PREPARE DATA
#LOAD MCX20
#LOAD MCX_COMBINED
#MCX20 METATABLE
#MCX_COMBINED METATABLE
#MCX20 DATATABLE
#MCX_COMBINED DATATABLE
#MCX20 DEG MARKERS
#MCX20 MSN MARKERS
#MCX20 FOXP2 MARKERS

#----FIG 4
#MCX_COMBINED CLUSTERPLOT
#MCX_COMBINED CLUSTERPLOT SPLIT BY DATASET
#PERCENTAGE BY DATASET
#HIERARCHICAL CLUSTERING
#DOTPLOT
#PCT MCX20

#----FIG 5
#MSN FOXP2
#MSN DR
#DR FOXP2
#MCX20_C_MSN_DR_FOXP2_pct_stackedbar
#DEG CORR PLOT

#----FIG 6
#CORRELATION PLOT
#RAINCLOUD PLOTS

#----SFIG 7
#TOPPGENE
#SFARI
#FOXP2 TARGETS

#----SFIG 8
#MCX_COMBINED DOTPLOT
#MCX20 CLUSTERPLOT
#MCX20 QC 
#MCX21 QC

#----SFIG 9
#HEATMAP SUPP

#----SFIG 10
#VIOLIN SUPP
#PN BLEND SUPPLEMENT

#LOAD MCX20
load(
  file = file.path(
    "MCX20_C_seurat.RData"
  )
)

#LOAD MCX_COMBINED
load(
  file = file.path(
  "MCX_COMBINED_C_seurat.RData"
  )
)

#MCX20 METATABLE
MCX20_C_cells <-
  enframe(Cells(MCX20_C_seurat), name = NULL) %>%
  rename("Cell" = value)

MCX20_C_clusters <-
  enframe(Idents(MCX20_C_seurat), name = NULL) %>%
  rename("Cluster" = value) 

MCX20_C_UMAP <-
  as_tibble(MCX20_C_seurat@reductions[["umap"]]@cell.embeddings)

MCX20_C_metadata <-
  FetchData(
    MCX20_C_seurat,
    vars = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mito",
      "orig.ident"
    )
  ) %>%
  as_tibble()

MCX20_C_metatable <-
  bind_cols(
    MCX20_C_cells,
    MCX20_C_metadata,
    MCX20_C_clusters,
    MCX20_C_UMAP
  ) %>%
  rename(
    UMIs = `nCount_RNA`,
    Genes = `nFeature_RNA`,
    pct_mito = percent.mito,
    Dataset = `orig.ident` 
  ) %>%
  add_count(Cluster, name = "n")

#MCX_COMBINED METATABLE
MCX_COMBINED_C_cells <-
  enframe(Cells(MCX_COMBINED_C_seurat), name = NULL) %>%
  rename("Cell" = value)

MCX_COMBINED_C_clusters <-
  enframe(Idents(MCX_COMBINED_C_seurat), name = NULL) %>%
  rename("Cluster" = value) 

MCX_COMBINED_C_UMAP <-
  as_tibble(MCX_COMBINED_C_seurat@reductions[["umap"]]@cell.embeddings)

MCX_COMBINED_C_metadata <-
  FetchData(
    MCX_COMBINED_C_seurat,
    vars = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mito",
      "orig.ident"
    )
  ) %>%
  as_tibble()

MCX_COMBINED_C_metatable <-
  bind_cols(
    MCX_COMBINED_C_cells,
    MCX_COMBINED_C_metadata,
    MCX_COMBINED_C_clusters,
    MCX_COMBINED_C_UMAP
  ) %>%
  rename(
    UMIs = `nCount_RNA`,
    Genes = `nFeature_RNA`,
    pct_mito = percent.mito,
    Dataset = `orig.ident` 
  ) %>%
  add_count(Cluster, name = "n")

#MCX20 DATATABLE
MCX20_C_data_tibble <-
  t(
    as.matrix(
      GetAssayData(MCX20_C_seurat, slot = "data", assay = "RNA")
    )
  ) %>%
  as_tibble(rownames = "Cell")

MCX20_C_cells <-
  enframe(Cells(MCX20_C_seurat), name = NULL) %>%
  rename("Cell" = value)

MCX20_C_clusters <-
  enframe(Idents(MCX20_C_seurat), name = NULL) %>%
  rename("Cluster" = value) 

MCX20_C_UMAP <-
  as_tibble(MCX20_C_seurat@reductions[["umap"]]@cell.embeddings)

MCX20_C_metadata <-
  FetchData(
    MCX20_C_seurat,
    vars = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mito",
      "seurat_clusters"
    )
  ) %>%
  as_tibble()

MCX20_C_datatable <-
  bind_cols(
    MCX20_C_cells,
    MCX20_C_metadata,
    MCX20_C_clusters,
    MCX20_C_UMAP
  ) %>%
  left_join(., MCX20_C_data_tibble, by = "Cell") %>%
  rename(
    UMIs = `nCount_RNA`,
    Genes = `nFeature_RNA`,
    `Percent Mito Genes` = percent.mito,
    `Cluster Number` = `seurat_clusters`
  ) %>%
  add_count(Cluster, name = "n") %>% 
  mutate(Cluster = as.numeric(as.character(Cluster)) + 1)

#MCX_COMBINED DATATABLE
MCX_COMBINED_C_data_tibble <-
  t(as.matrix(
    GetAssayData(MCX_COMBINED_C_seurat, slot = "data", assay = "RNA")
  )
  ) %>%
  as_tibble(rownames = "Cell")

MCX_COMBINED_C_cells <-
  enframe(Cells(MCX_COMBINED_C_seurat), name = NULL) %>%
  rename("Cell" = value)

MCX_COMBINED_C_clusters <-
  enframe(Idents(MCX_COMBINED_C_seurat), name = NULL) %>%
  rename("Cluster" = value) 

MCX_COMBINED_C_UMAP <-
  as_tibble(MCX_COMBINED_C_seurat@reductions[["umap"]]@cell.embeddings)

MCX_COMBINED_C_metadata <-
  FetchData(
    MCX_COMBINED_C_seurat,
    vars = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mito",
      "seurat_clusters",
      "orig.ident"
    )
  ) %>%
  as_tibble()

MCX_COMBINED_C_datatable <-
  bind_cols(
    MCX_COMBINED_C_cells,
    MCX_COMBINED_C_metadata,
    MCX_COMBINED_C_clusters,
    MCX_COMBINED_C_UMAP
  ) %>%
  left_join(., MCX_COMBINED_C_data_tibble, by = "Cell") %>%
  rename(
    UMIs = `nCount_RNA`,
    Genes = `nFeature_RNA`,
    `Percent Mito Genes` = percent.mito,
    `Cluster Number` = `seurat_clusters`,
    Dataset = `orig.ident` 
  ) %>%
  add_count(Cluster, name = "n")

#MCX20 DEG MARKERS
load(
  file = file.path(
    "MCX20_C_seurat_markers.RData"
  )
)

#MCX20 MSN MARKERS
load(
  file = file.path(
    "MCX20_C_seurat_MSN_markers.RData"
  )
)

#MCX20 FOXP2 MARKERS
load(
  file = file.path(
    "MCX_COMBINED_C_seurat_FOXP2_markers.RData"
  )
)

#MCX_COMBINED CLUSTERPLOT
MCX_COMBINED_C_seurat_clusterplot_labeled <-
  UMAPPlot(
    MCX_COMBINED_C_seurat,
    label = TRUE,
    label.size = 15,
    pt.size = 1.5,
    cols = rev(c(
      "red",
      "red",
      "red",
      "red",
      "red",
      "light blue",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "green",
      "green",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink"
    ))
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) 

ggsave(
  filename = file.path(
    "plot_clusterplot_combined.pdf"
  ),
  plot = MCX_COMBINED_C_seurat_clusterplot_labeled,
  width = 15,
  height = 15,
  units = "in",
  dpi = 300
) 

#MCX_COMBINED CLUSTERPLOT SPLIT BY DATASET
MCX_COMBINED_C_seurat_clusterplot_split <-
  UMAPPlot(
    MCX_COMBINED_C_seurat,
    split.by = "orig.ident",
    label = TRUE,
    label.size = 15,
    pt.size = 1.5,
    cols = rev(c(
      "red",
      "red",
      "red",
      "red",
      "red",
      "light blue",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "orange",
      "green",
      "green",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink"
    ))
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) 

ggsave(
  filename = file.path(
    "plot_clusterplot_split.pdf"
  ),
  plot = MCX_COMBINED_C_seurat_clusterplot_split,
  width = 30,
  height = 15,
  units = "in",
  dpi = 300
)

#PERCENTAGE BY DATASET
MCX_COMBINED_C_metatable_pct <-
  MCX_COMBINED_C_metatable %>%
  add_count(Cluster, Dataset, name = "nClusterDataset") %>%
  distinct(Cluster, Dataset, .keep_all = TRUE) %>% 
  mutate(pct = nClusterDataset / n * 100)

plot_MCX_COMBINED_C_metatable_pct <-
  MCX_COMBINED_C_metatable_pct %>%
  mutate(
    Cluster = factor(Cluster, levels = rev(c(
      2, 4, 1, 5, 3, 12, 13, 14, 16, 23, 19, 17, 
      18, 21, 9, 10, 6, 7, 8, 20, 22, 15, 11
    )))
  ) %>%
  ggplot(
    aes(
      x = Cluster,
      y = pct,
      fill = Dataset
    )
  ) +
  coord_flip() +
  geom_bar(
    stat = "identity", 
    position = position_stack(reverse = TRUE),
    color = "black"
  ) +
  scale_fill_manual(
    values = c("slategray1", "royalblue1"), 
    labels = c("shScr", "shFoxP2")
  ) +
  ylab("Percentage Composition") +
  scale_x_discrete(
    labels = rev(c(
      "2 (Medium Spiny Neuron)",
      "4 (Medium Spiny Neuron)",
      "1 (Medium Spiny Neuron)",
      "5 (Medium Spiny Neuron)",
      "3 (Medium Spiny Neuron)",
      "12 (Pallidal-like Neuron)",
      "13 (Interneuron)",
      "14 (Interneuron)",
      "16 (Interneuron)",
      "23 (Interneuron)",
      "19 (Interneuron)",
      "17 (Interneuron)",
      "18 (Interneuron)",
      "21 (Interneuron)",
      "9 (GLUT+ Neuron)",
      "10 (GLUT+ Neuron)",
      "6 (Astrocyte)",
      "7 (Astrocyte)",
      "8 (Oligodendrocyte)",
      "20 (Oligo. Precursor)",
      "22 (Microglia)",
      "15 (Endothelial)",
      "11 (Red Blood Cell)"
    ))
  ) +
  theme(
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    "plot_MCX_COMBINED_C_metatable_pct.pdf"
  ),
  plot = plot_MCX_COMBINED_C_metatable_pct,
  width = 9,
  height = 9,
  units = "in",
  dpi = 300
) 

#HIERARCHICAL CLUSTERING
MCX20_C_seurat_markers <-
  as_tibble(MCX20_C_seurat_markers) %>%
  mutate(cluster = as.numeric(cluster))

MCX20_C_seurat_avexp <- 
  log1p(AverageExpression(MCX20_C_seurat)$RNA) %>%
  as_tibble(rownames = "gene") 

nsig <- 
  50

markers_int_top <-
  MCX20_C_seurat_markers %>% 
  mutate(sign = avg_logFC > 0) %>%
  group_by(cluster, sign) %>%
  top_n(-1 * nsig, p_val_adj) %>%
  top_n(nsig, abs(avg_logFC)) %>%
  ungroup() %>% 
  distinct(gene, .keep_all = TRUE)

obj_int_avg_filt <-
  MCX20_C_seurat_avexp %>%
  filter(gene %in% markers_int_top$gene)

test_df <-
  as.data.frame(obj_int_avg_filt)

test_matrix <-
  data.matrix(test_df)[,-1]

rownames(test_matrix) <-
  test_df$gene

pv <-
  pvclust(
    test_matrix, 
    method.dist = "correlation",
    method.hclust = "ward.D",
    nboot = 200, 
    parallel = 4L
  )

plot(pv)

text_colors <-
  qualitative_hcl(6)

dend <-
  as.dendrogram(pv)

dend <-
  dend %>%
  dendextend::set("branches_k_color", value = text_colors, k = 6) 

dend_gg <-
  as.ggdend(dend)

dend_seg <-
  dend_gg$segments

scale_factor = .05
buffer_size = .1
position_box_width = 1
dend_seg = dend_seg %>%
  mutate(y2 = max(y)-y,
         yend2 = max(yend) - yend) %>%
  mutate(y2 = y2 * scale_factor,
         yend2 = yend2 * scale_factor) %>%
  mutate(col = if_else(is.na(col), "grey50", col))

tree_width = max(dend_seg$yend2)

dend_leaves <- 
  dend_gg$labels %>%
  mutate(
    test = case_when(
      label == 1 ~ "1 (Medium Spiny Neuron)",
      label == 5 ~ "5 (Medium Spiny Neuron)",
      label == 2 ~ "2 (Medium Spiny Neuron)",
      label == 3 ~ "3 (Medium Spiny Neuron)",
      label == 4 ~ "4 (Medium Spiny Neuron)",
      label == 10 ~ "10 (Pallidal Neuron)",
      label == 14 ~ "14 (Interneuron)",
      label == 15 ~ "15 (Interneuron)",
      label == 13 ~ "13 (Interneuron)", 
      label == 23 ~ "23 (Interneuron)",
      label == 8 ~ "8 (GLUT+ Unknown)",
      label == 11 ~ "11 (GLUT+ Interneuron)",
      label == 12 ~ "12 (Red Blood Cell)",
      label == 17 ~ "17 (Interneuron)",
      label == 19 ~ "19 (Interneuron)",
      label == 21 ~ "21 (Interneuron)",
      label == 22 ~ "22 (Interneuron)",
      label == 6 ~ "6 (Astrocyte)",
      label == 7 ~ "7 (Astrocyte)",
      label == 9 ~ "9 (Oligodendrocyte)",
      label == 20 ~ "20 (Oligodendrocyte Precursor)",
      label == 16 ~ "16 (Microglial)",
      label == 18 ~ "18 (Endothelial)"
    )
  )

dend_leaves = dend_leaves %>% 
  mutate(pos_xmax = x+.4,
         pos_xmin = x-.4,
         perc_ymin1 = min(dend_seg$yend) - buffer_size,
         perc_ymax1 = perc_ymin1 - position_box_width,
         perc_ymin2 = perc_ymax1 - buffer_size,
         perc_ymax2 = perc_ymin2 - position_box_width,
         label_y = perc_ymax2 - buffer_size 
  )

flat_plot <-
  ggplot() +
  geom_segment(data = dend_seg,
               aes(x = -x,
                   xend = -xend,
                   y = -y,
                   yend = -yend,
                   color = col
               ),
               lineend = "square",
               size = 1) + 
  coord_flip() +
  scale_color_identity() + 
  scale_fill_identity() +
  expand_limits(y=c(0, 5)) +
  theme_void() + 
  theme(legend.position = c(.1, .9)) +
  geom_text(data = dend_leaves,
            aes(x=-x,
                y=label_y,
                label=test),
            nudge_y = 2.5,
            hjust = 0
  )

ggsave(
  filename = file.path(
    "plot-tree-3.pdf"
  ),
  plot = flat_plot,
  width = 5,
  height = 8,
  units = "in",
  dpi = 300
)

#DOTPLOT
Idents(MCX20_C_seurat) <- 
  factor(
    Idents(MCX20_C_seurat),
    levels = rev(c(
      22, 19, 14, 17, 13, 21,
      23, 15, 11, 8, 10,
      7, 6, 
      1, 2, 3, 4, 5, 
      16, 12, 18, 9, 20
    ))
  )

dotplot <-
  DotPlot(
    MCX20_C_seurat,
    features = rev(c(
      "GAD2",
      "SLC17A6",
      "LHX6",
      "PVALB",
      "SST",
      "NPY",
      "NOS1",
      "CHAT",
      "PENK",
      "TSHZ1",
      "LRIG1",
      "MEIS2",
      "PPP1R1B",
      "FOXP1",
      "FOXP2",
      "TAC1",
      "CSF1R",
      "HBAA",
      "FLT1",
      "MBP",
      "PDGFRA"
    )),
    col.min = 0
  ) +
  scale_color_gradient(
    limits = c(0, 3),
    breaks = c(0, 3),
    low = "gray90",  
    high = "dodgerblue3"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, face = "italic"),
    legend.position = "top",
    legend.justification = "center"
  )

ggsave(
  filename = file.path(
    "dotplot.pdf"
  ),
  plot = dotplot,
  width = 5,
  height = 9.5,
  units = "in",
  dpi = 300
)

#PCT MCX20
MCX20_C_pct <-
  MCX20_C_metatable %>% 
  select(Cell, Cluster, n) %>%
  mutate(
    ClassSpecific = case_when(
      Cluster == 1 ~ "MSN-1",
      Cluster == 5 ~ "MSN-2",
      Cluster == 2 ~ "MSN-3",
      Cluster == 3 ~ "MSN-4",
      Cluster == 4 ~ "MSN-5",
      Cluster == 10 ~ "PN-1",
      Cluster == 14 ~ "INT-1",
      Cluster == 15 ~ "INT-2",
      Cluster == 13 ~ "INT-3",
      Cluster == 23 ~ "INT-4",
      Cluster == 8 ~ "GLUT-1",
      Cluster == 11 ~ "GLUT-2",
      Cluster == 12 ~ "RBC-1",
      Cluster == 17 ~ "INT-5",
      Cluster == 19 ~ "INT-6",
      Cluster == 21 ~ "INT-7",
      Cluster == 22 ~ "INT-8",
      Cluster == 6 ~ "A-1",
      Cluster == 7 ~ "A-2",
      Cluster == 9 ~ "O-1",
      Cluster == 20 ~ "OPC-1",
      Cluster == 16 ~ "M-1",
      Cluster == 18 ~ "E-1"
    ),
    ClassBroad = case_when(
      Cluster %in% c(1, 5, 2, 3, 4) ~ "MSN",
      Cluster %in% c(10) ~ "PN",
      Cluster %in% c(14, 15, 13, 23, 17, 19, 21, 22) ~ "INT",
      Cluster %in% c(8, 11) ~ "GLUT",
      Cluster %in% c(6, 7, 9, 20, 16, 18, 12) ~ "GLIAL"
    )
  ) %>%
  add_count(ClassSpecific, name = "ClassSpecific n") %>%
  add_count(ClassBroad, name = "ClassBroad n") %>%
  mutate(
    ClassSpecificPct = (`ClassSpecific n` / n()) * 100,
    ClassBroadPct = (`ClassBroad n` / n()) * 100
  ) %>%
  distinct(Cluster, .keep_all = TRUE) %>%
  arrange(ClassSpecific)

plot_pct <-
  MCX20_C_pct %>%
  distinct(ClassBroad, .keep_all = TRUE) %>%
  ggplot(
    aes(
      x = reorder(ClassBroad, ClassBroadPct),
      y = ClassBroadPct,
      fill = ClassBroad
    )
  ) +
  coord_flip() +
  geom_col() +
  scale_y_continuous(lim = c(0, 100)) +
  scale_fill_manual(
    values = c(
      "pink",
      "green",
      "orange",
      "red",
      "light blue"
    )
  ) +
  geom_text(
    aes(
      x = ClassBroad,
      y = ClassBroadPct, 
      label = round(ClassBroadPct, 2)
    ),
    hjust = 0,
    size = 24 / (14/5)
  ) +
  ylab("Percentage Composition") +
  scale_x_discrete(
    labels = rev(c(
      "Medium Spiny Neuron",
      "Glial Cell",
      "Interneuron",
      "GLUT+ Neuron", 
      "Pallidal-like Neuron"
    ))
  ) +
  theme(
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "plot_pct.pdf"
  ),
  plot = plot_pct,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

#MSN FOXP2
MCX20_C_FOXP2 <-
  MCX20_C_datatable %>%
  select(
    Cell,
    Cluster,
    UMAP_1,
    UMAP_2,
    FOXP2
  ) %>%
  mutate(
    Class = case_when(
      FOXP2 > 0 ~ "FOXP2+", FOXP2 == 0 ~ "FOXP2-"
    )
  ) 

MCX20_C_FOXP2 %>% 
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  mutate(nTotal = n()) %>%
  group_by(Class) %>% 
  mutate(
    Pct = round(n() / nTotal * 100, 0)
  ) %>%
  distinct(Class, Pct)

plot_MCX20_C_FOXP2 <-
  MCX20_C_FOXP2 %>%
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  ggplot(
    aes(x = UMAP_1, y = UMAP_2, color = Class)
  ) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(
      "orange", 
      "gray"
    ),
    limits = c(
      "FOXP2+",
      "FOXP2-"
    ),
    labels = c(
      "FoxP2+ (49%)",
      "FoxP2- (51%)"
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.025, 0.2),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_FOXP2.pdf"
  ),
  plot = plot_MCX20_C_FOXP2,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#MSN DR
MCX20_C_DR <-
  MCX20_C_datatable %>%
  mutate(
    Class = case_when(
      ((DRD1 > 0 | DRD5 > 0) & DRD2 > 0) ~ "DRD1/5+ & DRD2+",
      (DRD1 > 0 | DRD5 > 0) ~ "DRD1/5+",
      DRD2 > 0 ~ "DRD2+",
      (DRD1 == 0 & DRD5 == 0 & DRD2 == 0) ~ "DRD1- & DRD5- & DRD2-"
    )
  )


MCX20_C_DR %>% 
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  mutate(nTotal = n()) %>%
  group_by(Class) %>% 
  mutate(
    Pct = round(n() / nTotal * 100, 0)
  ) %>%
  distinct(Class, Pct)

plot_MCX20_C_DR <-
  MCX20_C_DR %>%
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  ggplot(
    aes(x = UMAP_1, y = UMAP_2, color = Class)
  ) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c("blue", "red", "green", "gray"), 
    limits = c("DRD1/5+", "DRD2+", "DRD1/5+ & DRD2+", "DRD1- & DRD5- & DRD2-"),
    labels = c("Drd1/5+ (36%)", "Drd2+ (14%)", "Drd1/5+ & Drd2+ (18%)", "Drd1- & Drd2- & Drd5- (33%)")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.02, 0.1),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_DR.pdf"
  ),
  plot = plot_MCX20_C_DR,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#DR FOXP2
MCX20_C_DR_FOXP2 <-
  MCX20_C_datatable %>%
  select(
    Cell,
    Cluster,
    UMAP_1,
    UMAP_2,
    DRD1,
    DRD5,
    DRD2,
    FOXP2
  ) %>%
  mutate(
    Class = case_when(
      (DRD1 > 0 | DRD5 > 0 | DRD2 > 0) &
        FOXP2 > 0 ~ "DR+ & FOXP2+",
      (DRD1 > 0 | DRD5 > 0) ~ "DRD1/5+",
      DRD2 > 0 ~ "DRD2+",
      FOXP2 > 0 ~ "FOXP2+",
      (DRD1 == 0 & DRD5 == 0 & DRD2 == 0) &
        FOXP2 == 0 ~ "DR- & FOXP2-"
    )
  ) 

MCX20_C_DR_FOXP2 %>% 
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  mutate(nTotal = n()) %>%
  group_by(Class) %>% 
  mutate(
    Pct = round(n() / nTotal * 100, 0)
  ) %>%
  distinct(Class, Pct)

plot_MCX20_C_DR_FOXP2 <-
  MCX20_C_DR_FOXP2 %>%
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  ggplot(
    aes(x = UMAP_1, y = UMAP_2, color = Class)
  ) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(
      "blue", 
      "red", 
      "orange",
      "green",
      "gray"
    ),
    limits = c(
      "DRD1/5+",
      "DRD2+",
      "FOXP2+",
      "DR+ & FOXP2+",
      "DR- & FOXP2-"
    ),
    labels = c(
      "Drd1/5+ (26%)", 
      "Drd2+ (11%)", 
      "FoxP2+ (19%)",
      "DR+ & FoxP2 (30%)",
      "DR- & FoxP2- (14%)"
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.025, 0.2),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_DR_FOXP2.pdf"
  ),
  plot = plot_MCX20_C_DR_FOXP2,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#MCX20_C_MSN_DR_FOXP2_pct_stackedbar
MCX20_C_MSN_DR_FOXP2_pct <-
  MCX20_C_DR %>% 
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -5 &
      `UMAP_2` > -6
  ) %>%
  mutate(
    Class = factor(
      Class, levels = c(
        "DRD1/5+",
        "DRD2+",
        "DRD1/5+ & DRD2+",
        "DRD1- & DRD5- & DRD2-"
      )
    ) 
  ) %>%
  mutate(nTotal = n()) %>%
  group_by(Class) %>%
  add_count(Class, `FOXP2+` = sum(FOXP2 > 0)) %>%
  mutate(
    Pct = n() / nTotal * 100,
    `PctFOXP2+` = `FOXP2+` / nTotal * 100
  ) %>%
  arrange(Class) %>%
  distinct(Class, .keep_all = TRUE) %>%
  pivot_longer(
    cols = c("Pct", "PctFOXP2+"),
    names_to = "type",
    values_to = "percentage"
  )

plot_MCX20_C_MSN_DR_pct_stackedbar <-
  MCX20_C_MSN_DR_FOXP2_pct %>%
  mutate(
    pctClass = round(`FOXP2+` / n * 100, 0),
    empty = 100 - pctClass
  ) %>%
  filter(type == "PctFOXP2+") %>%
  pivot_longer(
    cols = c("pctClass", "empty"),
    names_to = "group",
    values_to = "pct"
  ) %>% select(Class, group, pct) %>%
  ggplot(
    aes(
      x = fct_rev(Class),
      y = pct,
      fill = group
    )
  ) +
  geom_bar(stat = "identity", width = 0.25) +
  coord_flip() +
  labs(
    y = "Percentage of Cells"
  ) +
  scale_x_discrete(
    labels = rev(c(
      "Drd1/5+",
      "Drd2+",
      "Drd1/5+ & Drd2+",
      "Drd1/5- & Drd2-"
    ))
  ) +
  scale_fill_manual(
    values = rev(c(
      "orange1", 
      "gray40"
    )),
    labels = rev(c(
      "FoxP2+",
      "FoxP2-"
    )),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme(
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_MSN_DR_pct_stackedbar.pdf"
  ),
  plot = plot_MCX20_C_MSN_DR_pct_stackedbar,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

#DEG CORR PLOT

#MSN
MCX20_C_seurat_MSN <-
  subset(
    MCX20_C_seurat,
    idents = c(
      1, 5, 2, 3, 4
    ),
    UMAP_1 > -5 &
      UMAP_2 > -6
  )

MCX20_C_seurat_MSN_Group <-
  enframe(Idents(MCX20_C_seurat_MSN), name = NULL, value = "Cluster") %>%
  mutate(Group = case_when(
    Cluster %in% c(1, 5) ~ "A",
    Cluster %in% c(2, 3, 4) ~ "B"
  )
  ) %>%
  pull(Group)

MCX20_C_seurat_MSN$Group <- 
  MCX20_C_seurat_MSN_Group

Idents(MCX20_C_seurat_MSN) <- 
  "Group"

MCX20_C_seurat_MSN_avexp <- 
  log1p(
    AverageExpression(
      MCX20_C_seurat_MSN
    )
    $RNA
  ) %>%
  as_tibble(rownames = "gene") 

MCX20_C_seurat_MSN_markers_tibble <-
  as_tibble(MCX20_C_seurat_MSN_markers, rownames = "gene")

MCX20_C_seurat_MSN_avexp_stats <-
  inner_join(
    MCX20_C_seurat_MSN_avexp,
    MCX20_C_seurat_MSN_markers_tibble, 
    by = "gene"
  )

plot_MCX20_C_seurat_MSN_avexp_stats <-
  MCX20_C_seurat_MSN_avexp_stats %>%
  filter(
    gene %in% c(
      "ADARB2",
      "CHRM4",
      "DRD1",
      "DSCAM",
      "EBF1",
      "FOXP2",
      "ISL1",
      "LINGO2",
      "NRG3",
      "PDYN",
      "TAC1",
      "ADORA2A",
      "CHRM3",
      "DRD2",
      "GPR52",
      "GPR6",
      "PENK",
      "PTPRM"
    )
  ) %>%
  mutate(
    label = case_when(
      gene == "ADARB2" ~ "Adarb2",
      gene == "CHRM4" ~ "Chrm4",
      gene == "DRD1" ~ "Drd1",
      gene == "DSCAM" ~ "Dscam",
      gene == "EBF1" ~ "Ebf1",
      gene == "FOXP2" ~ "Foxp2",
      gene == "ISL1" ~ "Isl1",
      gene == "LINGO2" ~ "Lingo2",
      gene == "NRG3" ~ "Nrg3",
      gene == "PDYN" ~ "Pdyn",
      gene == "TAC1" ~ "Tac1",
      gene == "ADORA2A" ~ "Adora2a",
      gene == "CHRM3" ~ "Chrm3",
      gene == "DRD2" ~ "Drd2",
      gene == "GPR52" ~ "Gpr52",
      gene == "GPR6" ~ "Gpr6",
      gene == "PENK" ~ "Penk",
      gene == "PTPRM" ~ "Ptprm"
    ),
    highlight = case_when(
      gene %in% c(
        "ADARB2",
        "CHRM4",
        "DRD1",
        "DSCAM",
        "EBF1",
        "FOXP2",
        "ISL1",
        "LINGO2",
        "NRG3",
        "PDYN",
        "TAC1",
        "ADORA2A",
        "CHRM3",
        "DRD2",
        "GPR52",
        "GPR6",
        "PENK",
        "PTPRM"
      ) & p_val_adj >= 0.001 ~ "nonsig",
      gene %in% c(
        "ADARB2",
        "CHRM4",
        "DRD1",
        "DSCAM",
        "EBF1",
        "FOXP2",
        "ISL1",
        "LINGO2",
        "NRG3",
        "PDYN",
        "TAC1"
      ) & p_val_adj < 0.001 ~ "direct sig",
      gene %in% c(
        "ADORA2A",
        "CHRM3",
        "DRD2",
        "GPR52",
        "GPR6",
        "PENK",
        "PTPRM"
      ) & p_val_adj < 0.001 ~ "indirect sig"
    ),
    highlight = factor(
      highlight,
      levels = c(
        "nonsig",
        "direct sig",
        "indirect sig"
      )
    )
  ) %>%
  ggplot(
    aes(x = B, y = A, color = highlight)
  ) +
  geom_point(size = 3) +
  geom_abline(
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Indirect-like Pathway Average Expression",
    y = "Direct-like Pathway Average Expression"
  ) +
  geom_text_repel(
    aes(label = label, fontface = "bold.italic"), 
    size = (24 / (14/5)), 
    nudge_x = 0.75, 
    segment.size = 1.5, 
    show.legend = FALSE
  ) +
  lims(x = c(0, 4), y = c(0, 4)) +
  scale_color_manual(
    values = c("gray", "blue", "red"),
    labels = c("p > 0.001", "p < 0.001", "p < 0.001")
  ) +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.title.y = element_text(vjust = 2),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.position = c(0.75, 0.2)
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_seurat_MSN_avexp_stats.pdf"
  ),
  plot = plot_MCX20_C_seurat_MSN_avexp_stats,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#CORRELATION PLOT
DefaultAssay(MCX_COMBINED_C_seurat) <-
  "RNA"

MCX_COMBINED_C_seurat_MSN <-
  subset(
    MCX_COMBINED_C_seurat,
    idents = c(
      2, 4, 1, 5, 3
    ),
    UMAP_1 > -2 &
      UMAP_2 > -8
  )

DefaultAssay(MCX_COMBINED_C_seurat_MSN) <-
  "RNA"

MCX_COMBINED_C_seurat_MSN_FOXP2_DR <-
  subset(
    MCX_COMBINED_C_seurat_MSN,
    FOXP2 > 0 & (DRD1 > 0 | DRD2 > 0)
  )

Idents(MCX_COMBINED_C_seurat_MSN_FOXP2_DR) <- 
  MCX_COMBINED_C_seurat_MSN_FOXP2_DR$orig.ident

MCX_COMBINED_C_seurat_MSN_FOXP2_DR$D1 <- 
  FetchData(
    MCX_COMBINED_C_seurat_MSN_FOXP2_DR, 
    vars = c("DRD1")
  )

MCX_COMBINED_C_seurat_MSN_FOXP2_DR$D2 <- 
  FetchData(
    MCX_COMBINED_C_seurat_MSN_FOXP2_DR, 
    vars = c("DRD2")
  )

MCX_COMBINED_C_MSN_FOXP2_DR_D1D2corr <-
  as_tibble(
    list(
      DRD1 = MCX_COMBINED_C_seurat_MSN_FOXP2_DR$D1,
      DRD2 = MCX_COMBINED_C_seurat_MSN_FOXP2_DR$D2, 
      Dataset = MCX_COMBINED_C_seurat_MSN_FOXP2_DR$orig.ident
    )
  )

plot_MCX_COMBINED_C_MSN_FOXP2_DR_D1D2corr <- 
  MCX_COMBINED_C_MSN_FOXP2_DR_D1D2corr %>%
  ggplot(aes(
    x = DRD1, y = DRD2, color = Dataset
  )) +
  geom_jitter(
    width = 0.05,
    height = 0.05,
    size = 3
  ) +
  geom_abline(
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Drd1 Expression", 
    y = "Drd2 Expression"
  ) +
  scale_color_manual(
    labels = c("shScr", "shFoxP2"), 
    values = c("slategray2", "royalblue1")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 10)
    )
  ) +
  theme(
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 2),
    legend.text = element_text(size = 30),
    legend.position = "top",
    legend.title = element_blank(),
    legend.justification = "center",
    legend.spacing.x = unit(0.5, "cm")
  )

ggsave(
  filename = file.path(
    "plot_corr_test.pdf"
  ),
  plot = plot_MCX_COMBINED_C_MSN_FOXP2_DR_D1D2corr,
  width = 12,
  height = 12,
  units = "in",
  dpi = 300
)

#RAINCLOUD PLOTS
#https://github.com/RainCloudPlots/RainCloudPlots/blob/master/tutorial_R/R_rainclouds.R
#Allen M, Poggiali D, Whitaker K et al. Raincloud plots: a multi-platform tool for robust data visualization [version 2; peer review: 2 approved]. Wellcome Open Res 2021, 4:63. DOI: 10.12688/wellcomeopenres.15191.2

# Defining the geom_flat_violin function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )

MCX_COMBINED_C_MSN_test <-
  MCX_COMBINED_C_datatable %>%
  select(
    Cell,
    Cluster,
    Dataset,
    UMAP_1,
    UMAP_2,
    DRD1,
    DRD2,
    DRD5,
    FOXP2
  ) %>%
  filter(
    Cluster %in% c(
      1, 5, 2, 3, 4
    ) &
      `UMAP_1` > -2 &
      `UMAP_2` > -8 &
      FOXP2 > 0
  ) %>%
  mutate(
    Class = case_when(
      (DRD1 > 0 & DRD5 > 0 & DRD2 > 0) ~ "D1 & D5 & D2",
      (DRD1 > 0 & DRD5 > 0 & DRD2 == 0) ~ "D1 & D5",
      (DRD1 > 0 & DRD5 == 0 & DRD2 > 0) ~ "D1 & D2",
      (DRD1 == 0 & DRD5 > 0 & DRD2 > 0) ~ "D5 & D2",
      (DRD1 > 0 & DRD5 == 0 & DRD2 == 0) ~ "D1",
      (DRD1 == 0 & DRD5 == 0 & DRD2 > 0) ~ "D2",
      (DRD1 == 0 & DRD5 > 0 & DRD2 == 0) ~ "D5",
      (DRD1 == 0 & DRD5 == 0 & DRD2 == 0) ~ "D-"
    )
  )

plot_MCX_COMBINED_C_MSN_test_DRD1_raincloud <-
  MCX_COMBINED_C_MSN_test %>%
  mutate(
    Class = case_when(
      Class %in% c("D1", "D1 & D5") ~ "Direct",
      Class %in% c("D1 & D2", "D1 & D5 & D2") ~ "Mixed"
    ),
    Dataset = factor(
      Dataset,
      levels = c("MCX21", "MCX20")
    ),
    Class = factor(
      Class,
      levels = c("Mixed", "Direct")
    )
  ) %>%
  filter(Class != "NA") %>%
  ggplot(
    aes(
      x = Class, 
      y = DRD1,
      fill = Dataset
    )
  ) +
  geom_flat_violin(
    position = position_nudge(x = .2, y = 0), 
    adjust = 2, aes(alpha = 0.1)
  ) +
  scale_fill_manual(values = c("royalblue1", "slategray1"), guide = "none") +
  scale_alpha(guide = "none") +
  geom_point(
    position = position_jitter(width = .15),
    size = .25, 
    aes(color = Dataset, alpha = 0.1)
  ) +
  scale_color_manual(values = c("royalblue1", "slategray1"), 
                     labels = c("ShFoxP2", "ShScr")) +
  coord_flip() +
  geom_boxplot(
    aes(
      x = Class, 
      y = DRD1, 
      fill = Dataset,
      alpha = 0.1
    ), 
    outlier.shape = NA,
    width = .1, 
    colour = "black"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5, alpha = 0.9),
      reverse = TRUE
    )
  ) +
  ylim(NA, 3.25) +
  xlab("Drd1") +
  scale_x_discrete(labels = c("D1+/D2+ or \n D1+/D5+/D2+", "D1+ or \n D1+/D5+")) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(face = "italic"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

plot_MCX_COMBINED_C_MSN_test_DRD5_raincloud <-
  MCX_COMBINED_C_MSN_test %>%
  mutate(
    Class = case_when(
      Class %in% c("D5", "D1 & D5") ~ "Direct",
      Class %in% c("D5 & D2", "D1 & D5 & D2") ~ "Mixed"
    ),
    Dataset = factor(
      Dataset,
      levels = c("MCX21", "MCX20")
    ),
    Class = factor(
      Class,
      levels = c("Mixed", "Direct")
    )
  ) %>%
  filter(Class != "NA") %>%
  ggplot(
    aes(
      x = Class, 
      y = DRD5,
      fill = Dataset
    )
  ) +
  geom_flat_violin(
    position = position_nudge(x = .2, y = 0), 
    adjust = 2,
    aes(alpha = 0.1)
  ) +
  geom_point(
    position = position_jitter(width = .15),
    size = .25, 
    aes(color = Dataset, alpha = 0.1)
  ) +
  coord_flip() +
  geom_boxplot(
    aes(
      x = Class, 
      y = DRD5, 
      fill = Dataset, 
      alpha = 0.1
    ), 
    outlier.shape = NA,
    width = .1, 
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("royalblue1", "slategray2"), 
    guide = "none"
  ) +
  scale_color_manual(
    values = c("royalblue1", "slategray2")
  ) +
  scale_x_discrete(labels = c("D5+/D2+ or \n D1+/D5+/D2+", "D5+ or \n D1+/D5+")) +
  ylim(NA, 3.25) +
  xlab("Drd5") +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(face = "italic"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

plot_MCX_COMBINED_C_MSN_test_DRD2_raincloud <-
  MCX_COMBINED_C_MSN_test %>%
  mutate(
    Class = case_when(
      Class %in% c("D2") ~ "Indirect",
      Class %in% c("D1 & D2", "D5 & D2", "D1 & D5 & D2") ~ "Mixed"
    ),
    Dataset = factor(
      Dataset,
      levels = c("MCX21", "MCX20")
    ),
    Class = factor(
      Class,
      levels = c("Mixed", "Indirect")
    )
  ) %>%
  filter(Class != "NA") %>%
  ggplot(
    aes(
      x = Class, 
      y = DRD2,
      fill = Dataset
    )
  ) +
  geom_flat_violin(
    position = position_nudge(x = .2, y = 0), 
    adjust = 2, aes(alpha = 0.1)
  ) +
  geom_point(
    position = position_jitter(width = .15),
    size = .25, 
    aes(color = Dataset, alpha = 0.1)
  ) +
  coord_flip() +
  geom_boxplot(
    aes(
      x = Class, 
      y = DRD2, 
      fill = Dataset,
      alpha = 0.1
    ), 
    outlier.shape = NA,
    width = .1, 
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("royalblue1", "slategray2"), 
    guide = "none"
  ) +
  scale_color_manual(
    values = c("royalblue1", "slategray2")
  ) +
  scale_x_discrete(labels = c("D1+/D2+ or \n D5+/D2+ or \n D1+/D5+/D2+", "D2+")) +
  ylim(NA, 3.25) +
  xlab("Drd2") +
  ylab("Normalized Expression") +
  guides(
    color = guide_legend(
      override.aes = list(size = 5),
      reverse = TRUE
    )
  ) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "italic"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

testing <-
  ggarrange(
    plot_MCX_COMBINED_C_MSN_test_DRD1_raincloud,
    plot_MCX_COMBINED_C_MSN_test_DRD5_raincloud,
    plot_MCX_COMBINED_C_MSN_test_DRD2_raincloud
  )

ggsave(
  filename = file.path(
    "testing.pdf"
  ),
  plot = testing,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300
)

#TOPPGENE

#TOPPGENE OUTPUT DECREASED
toppgene_output_decreased <-
  read_tsv(
    file = file.path(
      "toppgene_output_decreased.txt"
    )
  ) %>%
  rename(
    qFDRBH = "q-value FDR B&H",
    Query = "Hit Count in Query List",
    Genome = "Hit Count in Genome",
    List = "Hit in Query List"
  ) %>%
  select(
    "Category",
    "ID",
    "Name",
    "qFDRBH",
    "Query",
    "Genome",
    "List"
  )

#TOPPGENE OUTPUT increased
toppgene_output_increased <-
  read_tsv(
    file = file.path(
      "toppgene_output_increased.txt"
    )
  ) %>%
  rename(
    qFDRBH = "q-value FDR B&H",
    Query = "Hit Count in Query List",
    Genome = "Hit Count in Genome",
    List = "Hit in Query List"
  ) %>%
  select(
    "Category",
    "ID",
    "Name",
    "qFDRBH",
    "Query",
    "Genome",
    "List"
  )

#REVIGO OUTPUT DECREASED
revigo_output_decreased_bp <-
  read_csv(
    file = file.path(
      "revigo_output_decreased_bp.csv"
    )
  ) %>%
  mutate(category = "Biological Process")

revigo_output_decreased_cc <-
  read_csv(
    file = file.path(
      "revigo_output_decreased_cc.csv"
    )
  ) %>%
  mutate(category = "Cellular Component")

revigo_output_decreased_mf <-
  read_csv(
    file = file.path(
      "revigo_output_decreased_mf.csv"
    )
  ) %>%
  mutate(category = "Molecular Function") 

revigo_output_decreased <-
  bind_rows(
    revigo_output_decreased_bp,
    revigo_output_decreased_cc,
    revigo_output_decreased_mf
  ) 

#TOPPGENE DECREASED FILTERED BY REVIGO DECREASED
toppgene_output_filtered_decreased <-
  toppgene_output_decreased %>%
  filter(
    ID %in% c(revigo_output_decreased$term_ID)
  ) %>%
  mutate(
    Change = "decrease"
  )

#REVIGO OUTPUT increased
revigo_output_increased_bp <-
  read_csv(
    file = file.path(
      "revigo_output_increased_bp.csv"
    )
  ) %>%
  mutate(category = "Biological Process")

revigo_output_increased_cc <-
  read_csv(
    file = file.path(
      "revigo_output_increased_cc.csv"
    )
  ) %>%
  mutate(category = "Cellular Component")

revigo_output_increased_mf <-
  read_csv(
    file = file.path(
      "revigo_output_increased_mf.csv"
    )
  ) %>%
  mutate(category = "Molecular Function") 

revigo_output_increased <-
  bind_rows(
    revigo_output_increased_bp,
    revigo_output_increased_cc,
    revigo_output_increased_mf
  )

#TOPPGENE INCREASED FILTERED BY REVIGO INCREASED
toppgene_output_filtered_increased <-
  toppgene_output_increased %>%
  filter(
    ID %in% c(revigo_output_increased$term_ID)
  ) %>%
  mutate(
    Change = "increase"
  )

#PLOT TOPPGENE
toppgene_output_merged <-
  bind_rows(
    toppgene_output_filtered_decreased,
    toppgene_output_filtered_increased
  ) %>%
  mutate(
    Category = factor(
      Category, levels = c(
        "GO: Biological Process",
        "GO: Cellular Component", 
        "GO: Molecular Function"
      )
    ),
    Name = fct_recode(
      Name, 
      "phosphorus transferase activity*" = 
        "transferase activity, transferring phosphorus-containing groups"
    )
  ) 

plot_toppgene_decreased <-
  toppgene_output_merged %>%
  filter(
    Change == "decrease"
  ) %>%
  group_by(
    Category
  ) %>%
  top_n(n = 5, Query) %>%
  arrange(Category, -Query) %>%
  rowid_to_column %>%
  ggplot(
    aes(
      x = reorder(Name, -rowid),
      y = Query,
      color = Category,
      label = Category
    )
  ) +
  coord_flip() +
  geom_point(
    size = 6, 
    position = position_dodge(0.75)
  ) +
  geom_linerange(
    aes(
      x = Name,
      ymin = 0,
      ymax = Query,
      color = Category
    ),
    size = 1.5, 
    alpha = 0.5, 
    position = position_dodge(0.75)
  ) +
  ggtitle("Decreased Expression in shFoxP2") +
  ylab("Number of Genes") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 36),
    axis.text.x = element_text(size = 36),
    axis.title.x = element_text(size = 36),
    plot.title = element_text(size = 36),
    text = element_text(size = 36),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank()
  )

plot_toppgene_increased <-
  toppgene_output_merged %>%
  filter(
    Change == "increase"
  ) %>%
  group_by(
    Category
  ) %>%
  top_n(n = 5, Query) %>%
  arrange(Category, -Query) %>%
  rowid_to_column %>%
  ggplot(
    aes(
      x = reorder(Name, -rowid),
      y = Query,
      color = Category,
      label = Category
    )
  ) +
  coord_flip() +
  geom_point(
    size = 6, 
    position = position_dodge(0.75)
  ) +
  geom_linerange(
    aes(
      x = Name,
      ymin = 0,
      ymax = Query,
      color = Category
    ),
    size = 1.5, 
    alpha = 0.5, 
    position = position_dodge(0.75)
  ) +
  ggtitle("Increased Expression in shFoxP2") +
  ylab("Number of Genes") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 36),
    axis.text.x = element_text(size = 36),
    axis.title.x = element_text(size = 36),
    plot.title = element_text(size = 36),
    text = element_text(size = 36),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank()
  )

plot_toppgene <-
  ggpubr::ggarrange(
    plot_toppgene_decreased,
    plot_toppgene_increased,
    nrow = 1,
    common.legend = TRUE,
    legend = "top"
  )

ggsave(
  filename = file.path(
    "plot_toppgene.pdf"
  ),
  plot = plot_toppgene,
  width = 24,
  height = 12,
  units = "in",
  dpi = 300
)

#SFARI

#LOAD SFARI
MCX_COMBINED_C_seurat_FOXP2_markers_sfari <-
  read_csv(
    file.path(
      "MCX_COMBINED_C_seurat_FOXP2_markers_sfari.csv"
    )
  )

#PLOT SFARI 
plot_sfari <-
  MCX_COMBINED_C_seurat_FOXP2_markers_sfari %>%
  filter(
    p_val_adj < 0.05,
    `gene-score` == "1"
  ) %>%
  select(-c(`gene-name`, `genetic-category`)) %>%
  mutate(
    flippedFC = avg_logFC * -1,
    change = case_when(
      avg_logFC < 0 ~ "increase",
      avg_logFC > 0 ~ "decrease"
    )
  ) %>%
  ggplot(
    aes(
      x = reorder(gene, -flippedFC),
      y = flippedFC,
      color = change
    )
  ) +
  coord_flip() +
  geom_point(
    size = 6, 
    position = position_dodge(0.75)
  ) +
  geom_linerange(
    aes(
      x = gene,
      ymin = 0,
      ymax = flippedFC,
      color = change
    ),
    size = 1.5, 
    alpha = 0.5, 
    position = position_dodge(0.75)
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed"
  ) +
  scale_color_discrete(
    labels = c("Decreased in shFoxP2", "Increased in shFoxP2"),
    name = "Expression Change"
  ) +
  ylab("Log Fold-Change") +
  scale_x_discrete(
    labels = rev(c(
      "Ankrd11",
      "Ptchd1",
      "Pcdh19",
      "Pten",
      "Ank3",
      "Hnrnpu",
      "Atrx",
      "Myt1l",
      "Arid1b",
      "Ank2",
      "Cacna1c",
      "Nsd1",
      "Gigyf2",
      "Ube3a",
      "Asxl3",
      "Bckdk",
      "Pogz",
      "Fmr1",
      "Nf1",
      "Vps13b",
      "Mbd5",
      "Nexmif",
      "Setd5",
      "Tcf4",
      "Phf21a",
      "Cdkl5",
      "Dyrk1a",
      "Arx",
      "Hras",
      "Nrxn3",
      "Rai1",
      "Baz2b",
      "Rfx3",
      "Aff2",
      "Arhgef9",
      "Phip",
      "Nrxn1",
      "Reln",
      "Nbea",
      "Katnal2",
      "Ddx3x",
      "Shank2"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    legend.position = c(1, 0.9),
    legend.justification = "right",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "plot_sfari.pdf"
  ),
  plot = plot_sfari,
  width = 10,
  height = 16,
  units = "in",
  dpi = 300
)

#FOXP2 TARGETS
foxp2_targets <-
  read_csv(
    file = file.path(
      "foxp2-targets.csv"
    )
  )

list <-
  foxp2_targets %>%
  pull(Gene)

targets <-
  MCX_COMBINED_C_seurat_FOXP2_markers %>%
  filter(gene %in% list, p_val_adj < 0.05) %>%
  arrange(avg_logFC) %>% print(n = Inf)

plot_targets <-
  targets %>%
  select(gene, avg_logFC) %>%
  mutate(
    flippedFC = avg_logFC * -1,
    change = case_when(
      avg_logFC < 0 ~ "increase",
      avg_logFC > 0 ~ "decrease"
    )
  ) %>%
  ggplot(
    aes(
      x = reorder(gene, -flippedFC),
      y = flippedFC,
      color = change
    )
  ) +
  coord_flip() +
  geom_point(
    size = 6, 
    position = position_dodge(0.75)
  ) +
  geom_linerange(
    aes(
      x = gene,
      ymin = 0,
      ymax = flippedFC,
      color = change
    ),
    size = 1.5, 
    alpha = 0.5, 
    position = position_dodge(0.75)
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed"
  ) +
  scale_color_discrete(
    labels = c("Decreased in shFoxP2", "Increased in shFoxP2"),
    name = "Expression Change"
  ) +
  ylab("Log Fold-Change") +
  scale_x_discrete(
    labels = rev(c(
      "Mras",
      "Snap25",
      "Fth1",
      "Cck",
      "Pdxk",
      "Ggfra1",
      "Gabbr1",
      "Ncf2",
      "Apod",
      "Mep1b",
      "Pax3",
      "Pnoc",
      "Dpagt1",
      "Plaa",
      "Ank1",
      "Cdc42bpb",
      "Fbxo22",
      "Nrip1",
      "Ppp2r5c",
      "Ndufa2",
      "Ppif",
      "Pkia",
      "Slc25a3",
      "Top2b",
      "Ubqln1",
      "Smoc1",
      "Akap6",
      "Fubp1",
      "Lipg",
      "Rps6ka2",
      "Dpp6",
      "Disc1",
      "Gpld1",
      "Ncoa1",
      "Calcrl",
      "Vldlr",
      "Mcf2",
      "Mef2d",
      "Mapre3",
      "Tgfb2",
      "Ptprm"
    ))
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    legend.position = c(1.1, 0.9),
    legend.justification = "right",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "plot_targets.pdf"
  ),
  plot = plot_targets,
  width = 10,
  height = 16,
  units = "in",
  dpi = 300
)

#MCX_COMBINED DOTPLOT
dotplot <-
  DotPlot(
    MCX_COMBINED_C_seurat,
    assay = "RNA",
    features = rev(c(
      "GAD2",
      "MEIS2",
      "PPP1R1B",
      "FOXP1",
      "FOXP2",
      "TAC1",
      "PENK",
      "TSHZ1",
      "LHX6",
      "PVALB",
      "SST",
      "NPY",
      "NOS1",
      "CHAT",
      "SLC17A6",
      "LRIG1",
      "MBP",
      "PDGFRA",
      "CSF1R",
      "FLT1",
      "HBAA"
    )),
    col.min = 0
  ) +
  scale_y_discrete(labels = rev(c(
    "2 (Medium Spiny Neuron)",
    "4 (Medium Spiny Neuron)",
    "1 (Medium Spiny Neuron)",
    "5 (Medium Spiny Neuron)",
    "3 (Medium Spiny Neuron)",
    "12 (Pallidal-like Neuron)",
    "13 (Interneuron)",
    "14 (Interneuron)",
    "16 (Interneuron)",
    "23 (Interneuron)",
    "19 (Interneuron)",
    "17 (Interneuron)",
    "18 (Interneuron)",
    "21 (Interneuron)",
    "9 (GLUT+ Neuron)",
    "10 (GLUT+ Neuron)",
    "6 (Astrocyte)",
    "7 (Astrocyte)",
    "8 (Oligodendrocyte)",
    "20 (Oligo. Precursor)",
    "22 (Microglia)",
    "15 (Endothelial)",
    "11 (Red Blood Cell)"
  ))) +
  scale_color_gradient(
    limits = c(0, 3),
    breaks = c(0, 3),
    low = "gray90",  
    high = "dodgerblue3"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, face = "italic"),
    legend.position = "top",
    legend.justification = "center"
  )

ggsave(
  filename = file.path(
    "dotplot.pdf"
  ),
  plot = dotplot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

#MCX20 CLUSTERPLOT
Idents(MCX20_C_seurat) <- 
  factor(
    Idents(MCX20_C_seurat),
    levels = rev(c(
      1, 5, 2, 3, 4, 10,
      14, 15, 13, 23,
      8,
      11, 12, 17, 19,
      21, 22,
      6, 7, 9, 20, 16, 18
    ))
  )

MCX20_C_seurat_clusterplot <-
  UMAPPlot(
    MCX20_C_seurat,
    label = TRUE,
    label.size = 10,
    pt.size = 1.5,
    cols = rev(c(
      "red",
      "red",
      "red",
      "red",
      "red",
      "light blue",
      "orange",
      "orange",
      "orange",
      "orange",
      "green",
      "green",
      "turquoise",
      "turquoise",
      "turquoise",
      "turquoise",
      "turquoise",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink",
      "pink"
    ))) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 30, vjust = -8),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_seurat_clusterplot_labeled.pdf"
  ),
  plot = MCX20_C_seurat_clusterplot,
  width = 15,
  height = 15,
  units = "in",
  dpi = 300
)

#MCX20 QC
plot_MCX20_C_seurat_UMI <-
  VlnPlot(
    MCX20_C_seurat,
    features = "nCount_RNA",
    pt.size = 0
  ) +
  geom_jitter(
    size = 2,
    alpha = 0.25
  ) +
  geom_hline(
    yintercept = 10000,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  ylim(-0.5, 50000) +
  labs(
    x = "Cell",
    y = "Number of UMI"
  ) +
  scale_x_discrete(labels = c(" ")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX20_C_seurat_mito <-
  VlnPlot(
    MCX20_C_seurat,
    features = "percent.mito",
    pt.size = 0
  ) +
  geom_jitter(
    size = 2,
    alpha = 0.25
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 20, 30, 40, 50, 60)
  ) +
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  labs(
    x = "Cell",
    y = "Percentage of Mitochondrial Genes"
  ) +
  scale_x_discrete(labels = c(" ")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX20_C_seurat_ngene <-
  FeatureScatter(
    MCX20_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    cols = "gray",
    pt.size = 2
  ) +
  xlim(0, 50000) +
  geom_segment(
    x = 10000,
    y = 0,
    xend = 10000,
    yend = 8000,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  scale_x_continuous(
    breaks = c(0, 10000, 20000, 30000, 40000, 50000),
    labels = c(0, "10k", "20k", "30k", "40k", "50k")
  ) +
  labs(
    x = "Number of UMI",
    y = "Number of Genes"
  ) +
  theme(
    axis.text.x = element_text(size = 32),
    axis.title.x = element_text(
      vjust = -1,
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX20_C_seurat_ncount <-
  FeatureScatter(
    MCX20_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "percent.mito",
    cols = "gray",
    pt.size = 2
  ) +
  xlim(0, 50000) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 20, 30, 40, 50, 60)
  ) +
  geom_segment(
    x = 0,
    y = 5,
    xend = 10000,
    yend = 5,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  geom_segment(
    x = 10000,
    y = 0,
    xend = 10000,
    yend = 5,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  scale_x_continuous(
    breaks = c(0, 10000, 20000, 30000, 40000, 50000),
    labels = c(0, "10k", "20k", "30k", "40k", "50k")
  ) +
  labs(
    x = "Number of UMI",
    y = "Percentage Mitochondrial Genes"
  ) +
  theme(
    axis.text.x = element_text(size = 32),
    axis.title.x = element_text(
      vjust = -1,
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX20_C_seurat_QC1 <-
  ggarrange(
    plot_MCX20_C_seurat_UMI,
    plot_MCX20_C_seurat_mito,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_seurat_QC1.pdf"
  ),
  plot = plot_MCX20_C_seurat_QC1,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300
)

plot_MCX20_C_seurat_QC2 <-
  ggarrange(
    plot_MCX20_C_seurat_ngene,
    plot_MCX20_C_seurat_ncount,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_seurat_QC2.pdf"
  ),
  plot = plot_MCX20_C_seurat_QC2,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300
)

#MCX21 QC
plot_MCX21_C_seurat_UMI <-
  VlnPlot(
    MCX21_C_seurat,
    features = "nCount_RNA",
    pt.size = 0
  ) +
  geom_jitter(
    size = 2,
    alpha = 0.25
  ) +
  geom_hline(
    yintercept = 10000,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  ylim(-0.5, 60000) +
  scale_y_continuous(
    breaks = c(0, 10000, 20000, 30000, 40000, 50000)
  ) +
  labs(
    x = "Cell",
    y = "Number of UMI"
  ) +
  scale_x_discrete(labels = c(" ")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX21_C_seurat_mito <-
  VlnPlot(
    MCX21_C_seurat,
    features = "percent.mito",
    pt.size = 0
  ) +
  geom_jitter(
    size = 2,
    alpha = 0.25
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 20)
  ) +
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  labs(
    x = "Cell",
    y = "Percentage of Mitochondrial Genes"
  ) +
  scale_x_discrete(labels = c(" ")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX21_C_seurat_ngene <-
  FeatureScatter(
    MCX21_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    cols = "gray",
    pt.size = 2
  ) +
  xlim(0, 60000) +
  scale_x_continuous(
    breaks = c(0, 10000, 20000, 30000, 40000, 50000) ,
    labels = c(0, "10k", "20k", "30k", "40k", "50k")
  ) +
  geom_segment(
    x = 10000,
    y = 0,
    xend = 10000,
    yend = 8000,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  labs(
    x = "Number of UMI",
    y = "Number of Genes"
  ) +
  theme(
    axis.text.x = element_text(size = 32),
    axis.title.x = element_text(
      vjust = -1,
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX21_C_seurat_ncount <-
  FeatureScatter(
    MCX21_C_seurat,
    feature1 = "nCount_RNA",
    feature2 = "percent.mito",
    cols = "gray",
    pt.size = 2
  ) +
  xlim(0, 60000) +
  scale_x_continuous(
    breaks = c(0, 10000, 20000, 30000, 40000, 50000),
    labels = c(0, "10k", "20k", "30k", "40k", "50k")
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 20)
  ) +
  geom_segment(
    x = 0,
    y = 5,
    xend = 10000,
    yend = 5,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  geom_segment(
    x = 10000,
    y = 0,
    xend = 10000,
    yend = 5,
    linetype = "dashed",
    size = 1,
    color = "red"
  ) +
  labs(
    x = "Number of UMI",
    y = "Percentage Mitochondrial Genes"
  ) +
  theme(
    axis.text.x = element_text(size = 32),
    axis.title.x = element_text(
      vjust = -1,
      size = 32
    ),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_text(
      vjust = 2,
      size = 32
    ),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot_MCX21_C_seurat_QC1 <-
  ggarrange(
    plot_MCX21_C_seurat_UMI,
    plot_MCX21_C_seurat_mito,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "plot_MCX21_C_seurat_QC1.pdf"
  ),
  plot = plot_MCX21_C_seurat_QC1,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300
)

plot_MCX21_C_seurat_QC2 <-
  ggarrange(
    plot_MCX21_C_seurat_ngene,
    plot_MCX21_C_seurat_ncount,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "plot_MCX21_C_seurat_QC2.pdf"
  ),
  plot = plot_MCX21_C_seurat_QC2,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300
)

#HEATMAP SUPP
MCX20_C_datatable <-
  MCX20_C_datatable %>%
  as_tibble() %>%
  mutate(
    Cluster = as.numeric(as.character(Cluster)) + 1,
    Cluster = factor(
      Cluster,
      levels = rev(c(
        1, 5, 2, 3, 4,
        10, 
        14, 15, 13, 23, 12, 17, 19, 21, 22,
        8, 11,
        6, 7, 9, 20, 16, 18
      ))
    )
  )

MCX20_C_heatmap_supp <-
  MCX20_C_datatable %>%
  select(
    Cluster,
    CBLN1,
    LHX1,
    PVALB,
    SST,
    GREM1,
    MEIS2,
    `NKX2-1`,
    PENK,
    PKIB,
    SCN4B,
    CALB2,
    DEPTOR,
    ELFN1,
    GRIK3,
    LHX6,
    PAPPA,
    VIP,
    FIBCD1,
    GBX1,
    NPR3,
    CASZ1,
    TRDN,
    FOXP1,
    FOXP2,
    SYT1,
    BARHL1,
    FOXA1,
    LMX1A,
    PITX2
  ) %>%
  pivot_longer(
    cols = -1,
    names_to = "Gene", 
    values_to = "Expression"
  ) %>%
  group_by(Gene, Cluster) %>%
  mutate(
    Expression = mean(Expression)
  ) %>%
  distinct(Cluster, .keep_all = TRUE)

plot_MCX20_C_heatmap_supp_CBLN1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "CBLN1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_x_discrete(labels = c(
    "1", "5", "2", "3", "4",
    "10", 
    "14", "15", "13", "23", "12", "17", "19", "21", "22",
    "8", "11",
    "6", "7", "9", "20", "16", "18"
  ), position = "top") +
  scale_y_discrete(labels = c("Cbln1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_LHX1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "LHX1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Lhx1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_PVALB <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "PVALB") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Pvalb")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_SST <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "SST") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Sst")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_GREM1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "GREM1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Grem1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_MEIS2 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "MEIS2") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  ) +
  scale_y_discrete(labels = c("Meis2")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_NKX2.1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "NKX2-1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Nkx2-1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_PENK <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "PENK") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Penk")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_PKIB <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "PKIB") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Pkib")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_SCN4B <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "SCN4B") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Scn4b")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_CALB2 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "CALB2") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Calb2")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_DEPTOR <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "DEPTOR") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Deptor")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_ELFN1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "ELFN1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_y_discrete(labels = c("Elfn1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_GRIK3 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "GRIK3") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Grik3")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_LHX6 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "LHX6") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Lhx6")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_PAPPA <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "PAPPA") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 2),
    breaks = c(0, 2)
  ) +
  scale_y_discrete(labels = c("Pappa")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_VIP <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "VIP") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Vip")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_FIBCD1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "FIBCD1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Fibcd1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_GBX1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "GBX1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Gbx1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_NPR3 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "NPR3") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Npr3")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_CASZ1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "CASZ1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Casz1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_TRDN <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "TRDN") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Trdn")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_FOXP1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "FOXP1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  ) +
  scale_y_discrete(labels = c("FoxP1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_FOXP2 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "FOXP2") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 3),
    breaks = c(0, 3)
  ) +
  scale_y_discrete(labels = c("FoxP2")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_SYT1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "SYT1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 4),
    breaks = c(0, 4)
  ) +
  scale_y_discrete(labels = c("Syt1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_BARHL1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "BARHL1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Barhl1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_FOXA1 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "FOXA1") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("FoxA1")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_LMX1A <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "LMX1A") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Lmx1a")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp_PITX2 <-
  MCX20_C_heatmap_supp %>%
  filter(Gene == "PITX2") %>%
  ggplot(aes(x = fct_rev(Cluster), y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient(
    low = "white",  
    high = "steelblue",
    limits = c(0, 1),
    breaks = c(0, 1)
  ) +
  scale_y_discrete(labels = c("Pitx2")) +
  labs(fill = element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 24),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

plot_MCX20_C_heatmap_supp <-
  ggarrange(
    plot_MCX20_C_heatmap_supp_CBLN1,
    plot_MCX20_C_heatmap_supp_LHX1,
    plot_MCX20_C_heatmap_supp_FOXP2,
    plot_MCX20_C_heatmap_supp_PVALB,
    plot_MCX20_C_heatmap_supp_SST,
    plot_MCX20_C_heatmap_supp_FOXP2,
    plot_MCX20_C_heatmap_supp_GREM1,
    plot_MCX20_C_heatmap_supp_MEIS2,
    plot_MCX20_C_heatmap_supp_NKX2.1,
    plot_MCX20_C_heatmap_supp_PENK,
    plot_MCX20_C_heatmap_supp_PKIB,
    plot_MCX20_C_heatmap_supp_PVALB,
    plot_MCX20_C_heatmap_supp_SCN4B,
    plot_MCX20_C_heatmap_supp_PITX2,
    plot_MCX20_C_heatmap_supp_CALB2,
    plot_MCX20_C_heatmap_supp_DEPTOR,
    plot_MCX20_C_heatmap_supp_ELFN1,
    plot_MCX20_C_heatmap_supp_GRIK3,
    plot_MCX20_C_heatmap_supp_LHX6,
    plot_MCX20_C_heatmap_supp_PAPPA,
    plot_MCX20_C_heatmap_supp_VIP,
    plot_MCX20_C_heatmap_supp_FIBCD1,
    plot_MCX20_C_heatmap_supp_GBX1,
    plot_MCX20_C_heatmap_supp_NPR3,
    plot_MCX20_C_heatmap_supp_CASZ1,
    plot_MCX20_C_heatmap_supp_TRDN,
    plot_MCX20_C_heatmap_supp_FOXP1,
    plot_MCX20_C_heatmap_supp_FOXP2,
    plot_MCX20_C_heatmap_supp_FOXA1,
    plot_MCX20_C_heatmap_supp_LMX1A,
    plot_MCX20_C_heatmap_supp_BARHL1,
    plot_MCX20_C_heatmap_supp_SYT1,
    ncol = 1
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_heatmap_supp.pdf"
  ),
  plot = plot_MCX20_C_heatmap_supp,
  width = 24,
  height = 29,
  units = "in",
  limitsize = FALSE,
  dpi = 300
)

#VIOLIN SUPP
pn_violin <- 
  c(
    "MEIS2",
    "FOXP2",
    "PENK",
    "PVALB",
    "NKX2-1",
    "LHX6",
    "GREM1",
    "SCN4B",
    "CBLN1",
    "LHX1",
    "SST"
  )

pn_violin_plots <- 
  lapply(
    pn_violin, function(x) {
      VlnPlot(
        MCX20_C_seurat_PN,
        features = x,
        pt.size = 0,
        assay = "RNA",
        slot = "data"
      )
    } +
      coord_flip() +
      theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 24, angle = 0),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold.italic"),
        legend.position = "none"
      )
  )

names(pn_violin_plots) <-
  pn_violin

plot_pn_violin_MEIS2 <-
  pn_violin_plots$MEIS2 +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  ggtitle("Meis2") +
  scale_fill_manual(values = "firebrick") +
  theme(
    axis.text.y = element_blank()
  )

plot_pn_violin_FOXP2 <-
  pn_violin_plots$FOXP2 +
  scale_y_continuous(
    limits = c(NA, 5),
    breaks = c(0, 5)
  ) +
  ggtitle("FoxP2") +
  scale_fill_manual(values = "firebrick") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_PENK <-
  pn_violin_plots$PENK +
  scale_y_continuous(
    limits = c(NA, 6),
    breaks = c(0, 6)
  ) +
  ggtitle("Penk") +
  scale_fill_manual(values = "firebrick") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_PVALB <-
  pn_violin_plots$PVALB +
  scale_y_continuous(
    limits = c(NA, 3),
    breaks = c(0, 3)
  ) +
  ggtitle("Pvalb") +
  scale_fill_manual(values = "green") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_NKX2.1 <-
  pn_violin_plots$"NKX2-1" +
  scale_y_continuous(
    limits = c(NA, 3),
    breaks = c(0, 3)
  ) +
  ggtitle("Nkx2-1") +
  scale_fill_manual(values = "green") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_LHX6 <-
  pn_violin_plots$LHX6 +
  scale_y_continuous(
    limits = c(NA, 3),
    breaks = c(0, 3)
  ) +
  ggtitle("Lhx6") +
  scale_fill_manual(values = "green") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_GREM1 <-
  pn_violin_plots$GREM1 +
  scale_y_continuous(
    limits = c(NA, 3),
    breaks = c(0, 3)
  ) +
  ggtitle("Grem1") +
  scale_fill_manual(values = "green") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_SCN4B <-
  pn_violin_plots$SCN4B +
  scale_y_continuous(
    limits = c(NA, 2),
    breaks = c(0, 2)
  ) +
  ggtitle("Scn4b") +
  scale_fill_manual(values = "green") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_CBLN1 <-
  pn_violin_plots$CBLN1 +
  scale_y_continuous(
    limits = c(NA, 2),
    breaks = c(0, 2)
  ) +
  ggtitle("Cbln1") +
  scale_fill_manual(values = "royal blue") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_LHX1 <-
  pn_violin_plots$LHX1 +
  scale_y_continuous(
    limits = c(NA, 1),
    breaks = c(0, 1)
  ) +
  ggtitle("Lhx1") +
  scale_fill_manual(values = "royal blue") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_PVALB2 <-
  pn_violin_plots$PVALB +
  scale_y_continuous(
    limits = c(NA, 3),
    breaks = c(0, 3)
  ) +
  ggtitle("Pvalb") +
  scale_fill_manual(values = "royal blue") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_pn_violin_SST <-
  pn_violin_plots$SST +
  scale_y_continuous(
    limits = c(NA, 3),
    breaks = c(0, 3)
  ) +
  ggtitle("Sst") +
  scale_fill_manual(values = "royal blue") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )

vlnplot <-
  ggarrange(
    plot_pn_violin_MEIS2,
    plot_pn_violin_FOXP2,
    plot_pn_violin_PENK,
    plot_pn_violin_PVALB,
    plot_pn_violin_NKX2.1,
    plot_pn_violin_LHX6,
    plot_pn_violin_GREM1,
    plot_pn_violin_SCN4B,
    plot_pn_violin_CBLN1,
    plot_pn_violin_LHX1,
    plot_pn_violin_PVALB2,
    plot_pn_violin_SST,
    nrow = 1
  )

ggsave(
  filename = file.path(
    "vlnplot.pdf"
  ),
  plot = vlnplot,
  width = 16,
  height = 4,
  units = "in",
  dpi = 300
)

#PN BLEND SUPPLEMENT
MCX20_C_seurat_PN <-
  subset(
    MCX20_C_seurat,
    idents = 10,
    UMAP_1 > 0 &
      UMAP_2 > -7.5
  )

MCX20_C_seurat_PN_blend <-
  FeaturePlot(
    MCX20_C_seurat_PN, 
    features = c("FOXP2", "PENK"), 
    blend = TRUE,
    combine = FALSE
  )

plot_MCX20_C_seurat_PN_blend <-
  MCX20_C_seurat_PN_blend[[3]] +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_seurat_PN_blend.pdf"
  ),
  plot = plot_MCX20_C_seurat_PN_blend,
  width = 8,
  height = 8,
  units = "in",
  limitsize = FALSE,
  dpi = 300
)

plot_MCX20_C_seurat_PN_blend_legend <-
  MCX20_C_seurat_PN_blend[[4]] +
  labs(
    x = "Relative FoxP2 Expression",
    y = "Relative Penk Expression"
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = c(0, 10)
  ) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = c(0, 10)
  ) +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24)
  )

ggsave(
  filename = file.path(
    "plot_MCX20_C_seurat_PN_blend_legend.pdf"
  ),
  plot = plot_MCX20_C_seurat_PN_blend_legend,
  width = 4,
  height = 4,
  units = "in",
  limitsize = FALSE,
  dpi = 300
)
