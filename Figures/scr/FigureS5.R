library(Seurat)
library(future)
plan("multisession", workers = 2)
library(hdf5r)
library(ggplot2)
library(spacexr)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(patchwork)
options(future.globals.maxSize = 2 * 1024^3)
#BiocManager::install("glmGamPoi")
library(glmGamPoi)
# Load packages
Sys.setenv(OPENAI_API_KEY = '######') # your API key
library(GPTCelltype)
library(openai)

options(bitmapType = "cairo")


setwd("/data/l2_bioinfo1/lyul1/IBD/sp_transcriptomic/spatial")
dge_analysis_plot<- function(sample) {
  
  spatialobj<-Load10X_Spatial(paste0("./", sample), slice="S1") # a Seurat object containing gene expression data and spatial coordinates.
  spatialobj<- spatialobj %>% 
    SCTransform(assay = "Spatial", verbose = FALSE)  %>% # Normalizes gene expression (log normalization).
    FindVariableFeatures() %>% # Identifies highly variable genes (most informative for downstream analysis).
    ScaleData() %>% # Scales and centers gene expression values.
    RunPCA() %>% # Performs Principal Component Analysis (PCA) for dimensionality reduction.
    RunUMAP(dims=1:30) %>% # Runs Uniform Manifold Approximation and Projection (UMAP) for visualization in lower dimensions.
    FindNeighbors() %>% # Computes a nearest-neighbor graph based on PCA.
    FindClusters() # Groups similar cells into clusters using Louvain algorithm.
  
  # plot
  # umap
  p <- DimPlot(spatialobj, reduction = "umap")
  ggsave(paste0("/data/l2_bioinfo1/lyul1/IBD/sp_transcriptomic/plots/umap_plot_", sample, ".pdf"), plot = p, width = 6, height = 4, dpi = 300)
  
  # Data preprocessing
  plot1 <- VlnPlot(spatialobj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(spatialobj, features = "nCount_Spatial") + theme(legend.position = "right")
  p <- wrap_plots(plot1, plot2)
  ggsave(paste0("/data/l2_bioinfo1/lyul1/IBD/sp_transcriptomic/plots/data_prepro_", sample, ".pdf"), plot = p, width = 6, height = 4, dpi = 300)
  
  # Dimensionality reduction, clustering, and visualization
  # We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.
  # We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().
  p1 <- DimPlot(spatialobj, reduction = "umap", label = TRUE, label.size = 3)
  p2 <- SpatialDimPlot(spatialobj, label = TRUE, label.size = 3)
  p3<- p1 + p2
  #SpatialDimPlot(spatialobj, cells.highlight = CellsByIdentities(object = spatialobj, idents = c(2, 1, 4, 3,
  #                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)
    ggsave(paste0("/data/l2_bioinfo1/lyul1/IBD/sp_transcriptomic/plots/SpatialDimPlotoverlaid_", sample, ".pdf"), plot = p3, width = 6, height = 4, dpi = 300)

### annotation using "GPTCelltype"  (Hou, W., Ji, Z. Assessing GPT-4 for cell type annotation in single-cell RNA-seq analysis. Nat Methods 21, 1462–1465 (2024).)

  ### GPT-based annotations may not be fully reproducible because API outputs can vary across runs and model updates, even with identical inputs.
  markers <- FindAllMarkers(spatialobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  res <- gptcelltype(markers, model = 'gpt-4', 
                     tissue = "colon")
  # Assign cell type annotation back to Seurat object
  spatialobj@meta.data$celltype <- as.factor(res[as.character(Idents(spatialobj))])
  write.csv(spatialobj@meta.data, paste0("/data/l2_bioinfo1/lyul1/IBD/sp_transcriptomic/results/", "gpt_cell_types_colon_pos_", sample, ".csv"), row.names = F)
