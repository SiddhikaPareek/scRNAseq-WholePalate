---
  title: "Pseudotime analysis"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggrepel)
library(data.table)
library(GeneSwitches)
library(SingleCellExperiment)

remotes::install_github('satijalab/seurat-wrappers')


Integration.mesenchyme <- subset(Integration, idents = c("3","0","2","1","4","12","5","9","7"))

cds <- as.cell_data_set(Integration.mesenchyme)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


FeaturePlot(object = Integration.mesenchyme, features = c("Hic1", "Piezo2","Runx2","Tbx22","Tfap2b","Meox2", "Shox2", "Top2a", "Dlx5", "Sp7", "Tubb3", "Wnt16", "Tbx22", "Msx1", "Sox9", "Pax9", "Wwp2", "Ugdh", "Susd5", "Sox10", ), min.cutoff = "q9", sort.cell = "true")
FeaturePlot(object = Integration.mesenchyme, features = c( "Hic1", "Piezo2"), min.cutoff = "q9", sort.cell = "true")

#### DEBUG PSEUDOTIME MONOCLE 3 FROM SEURAT METADATA ####

rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds@reduce_dim_aux$UMAP) <- NULL
colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL

cds.new <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_sub <-cds
cds_sub <- choose_graph_segments(cds,clear_cds = F)
cds_sub <- choose_graph_segments(cds_from_seurat,clear_cds = F)
cds_sub  <- preprocess_cds(cds_sub ,method = "PCA", num_dim = 50)
cds_sub  <- reduce_dimension(cds_sub)
cds_sub <- cluster_cells(cds = cds_sub, reduction_method = "UMAP")
cds_sub <- learn_graph(cds_sub)
cds_sub  <- order_cells(cds_sub)

plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds_sub, genes=c("Axin2", "Aox3", "Dspp", "Crabp1","Tac1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,cell_size = 0.5)
cds_sub2 <- choose_graph_segments(cds_from_seurat,clear_cds = F)
cds_sub  <- preprocess_cds(cds_sub ,method = "PCA", num_dim = 50)
cds_sub  <- reduce_dimension(cds_sub)
cds_sub <- cluster_cells(cds = cds_sub, reduction_method = "UMAP")
cds_sub <- learn_graph(cds_sub)
cds_sub2  <- order_cells(cds_sub2)

plot_cells(cds_sub2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

try1 <- CreateSeuratObject(
  cds_sub@assays@data@listData$counts,
  project = "CreateSeuratObject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)

# STEP 2
try1 <- NormalizeData(try1, normalization.method = "LogNormalize", scale.factor = 10000)

# STEP 3 + removing genes with zero variance
data <- as.matrix(try1@assays$RNA@data)
data_filtered <- data[apply(data[,-1], 1, function(x) !all(x==0)),]

##Isolate pseudotime and umap data
Pseudotime <- cds_sub@principal_graph_aux@listData$UMAP$pseudotime
UMAP <- cds_sub@int_colData@listData$reducedDims$UMAP
str(object = cds_sub)
sce_object <- SingleCellExperiment(assays = List(expdata = data_filtered))
colData(sce_object)$Pseudotime <- Pseudotime
reducedDims(sce_object) <- SimpleList(UMAP = UMAP)


#From here you can follow the GeneSwitch tutorial: 
# https://geneswitches.ddnetbio.com/

### check the threshold for fast binarization
h <- hist(assays(sce_object)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.2, col="blue")}


bn_cutoff <- 0.7
sce_object <- binarize_exp(sce_object, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff) 
sce_object <- find_switch_logistic_fastglm(sce_object, downsample = TRUE, show_warning = FALSE)

sg_allgenes <- filter_switchgenes(sce_object, allgenes = TRUE, topnum = 200) 

length(sg_allgenes@listData[["geneID"]])

## filter top 15 best fitting switching genes among all the genes
sg_allgenes <- filter_switchgenes(sce_object, allgenes = TRUE, topnum = 30)
## filter top 15 best fitting switching genes among surface proteins and TFs only
sg_gtypes <- filter_switchgenes(sce_object, allgenes = FALSE, topnum = 30,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))
## combine switching genes and remove duplicated genes from sg_allgenes
sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])

plot_timeline_ggplot(sg_vis, timedata = sce_object$Pseudotime ,txtsize = 5)

plot_gene_exp(sce_object, gene = "VIM", reduction = "monocleRD", downsample = F)
