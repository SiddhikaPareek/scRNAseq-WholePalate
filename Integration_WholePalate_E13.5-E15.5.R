# Assigning stage to Mesenchyme subcluster seurat object (to be used as input forintegration)
E13.5whole$stage <- "E13.5"
E14.5whole$stage <- "E14.5"
E15.5whole$stage <- "E15.5"

# Perform integration   
MesenchymeIntegrated.anchors <- FindIntegrationAnchors(object.list = list(E13.5whole, E14.5whole, E15.5whole), anchor.features = 5000, dims = 1:20)
MesenchymeIntegrated.combined <- IntegrateData(anchorset = MesenchymeIntegrated.anchors, dims = 1:20)
dim(MesenchymeIntegrated.combined) # cells=2000 genes=14478

DefaultAssay(object = MesenchymeIntegrated.combined) <- "integrated"

library(dplyr)
library(stringr)

# Regress cellcycle genes
mouse_cell_cycle_genes <- readRDS("C:/Users/siddh/OneDrive - University of Southern California/Palate/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes[[1]]
g2m.genes <- mouse_cell_cycle_genes[[2]]
MesenchymeIntegrated.combined <- CellCycleScoring(MesenchymeIntegrated.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
MesenchymeIntegrated.combined <- RunPCA(MesenchymeIntegrated.combined, features = c(s.genes, g2m.genes))
  # ^ this line is creating PCs using just the cell cycle genes 

p1 <- DimPlot(MesenchymeIntegrated.combined)
MesenchymeIntegrated.combined <- ScaleData(MesenchymeIntegrated.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(MesenchymeIntegrated.combined))

MesenchymeIntegrated.combined <- FindVariableFeatures(MesenchymeIntegrated.combined, selection.method = "vst", nfeatures = 500)
# ^ this will find new variables features for the integrated object 

MesenchymeIntegrated.combined <- RunPCA(MesenchymeIntegrated.combined, npcs = 50)
#MesenchymeIntegrated.combined <- RunPCA(MesenchymeIntegrated.combined, features = c(s.genes, g2m.genes))

ElbowPlot(MesenchymeIntegrated.combined, ndims = 50)


p2 <- DimPlot(MesenchymeIntegrated.combined)

p1+p2
plot(p2)

# UMAP and Clustering for integrted object
MesenchymeIntegrated.combined <- RunUMAP(object = MesenchymeIntegrated.combined, reduction = "pca", dims = 1:10)
  # ^ update dims based on elbow plot

MesenchymeIntegrated.combined <- FindNeighbors(object = MesenchymeIntegrated.combined, reduction = "pca", dims = 1:10)
  # ^ update dims based on elbow plot

MesenchymeIntegrated.combined <- FindClusters(MesenchymeIntegrated.combined, resolution = 0.5)
  # ^ try various resolutions 0.5, 0.7, 1.0 etc. and see how the clusters appear on the UMAPs 

saveRDS(MesenchymeIntegrated.combined,file="C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/Integrated/MesenchymeIntegrated.combined.rds")

# Visualization
p1 <- DimPlot(object = MesenchymeIntegrated.combined, reduction = "umap", label = TRUE, split.by = "stage") +
  ggtitle("All Cells - Colored by Seurat4 Clustering") + theme(plot.title = element_text(hjust=0.5)) +
  xlim(-12,13) + ylim(-14,17)
p2 <- DimPlot(object = MesenchymeIntegrated.combined, reduction = "umap", group.by = "stage") +
  ggtitle("All Cells - Colored by Cell Stages") + theme(plot.title = element_text(hjust=0.5)) +
  xlim(-12,10) + ylim(-16,13)
p4 <- DimPlot(object = MesenchymeIntegrated.combined, reduction = "umap", group.by = "Phase")+ 
  ggtitle("All Cells - Colored by Seurat3 \nCell Cycle Phase") + theme(plot.title = element_text(hjust=0.5)) +
  xlim(-12,10) + ylim(-16,13)
plot(p1)
plot(p2)
plot(p4)
par(mfrow = c(1, 1),oma=c(1,1,1,1),mar=c(4,4,4,4)+0.1) 
graphics.off()

# Identify conserved cell type markers
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')

#DefaultAssay(object = MesenchymeIntegrated.combined) <- "RNA"
#nk.markers <- FindConservedMarkers(object = MesenchymeIntegrated.combined, ident.1 = 11, grouping.var = "stage", 
#                                   verbose = FALSE)
#head(x = nk.markers)
#write.csv(nk.markers, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/Integrated/Integrated_ConservedMarkers.csv")
#FeaturePlot(object = MesenchymeIntegrated.combined, features = c("Hic1", "Piezo2","Runx2","Tbx22","Tfap2b","Meox2", "Shox2", "Top2a", "Dlx5", "Sp7", "Tubb3", "Wnt16", "Tbx22", "Msx1", "Sox9", "Pax9", "Wwp2", "Ugdh", "Susd5", "Sox10"), min.cutoff = "q9")

# Identify conserved cell type markers
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
library(metap)
DefaultAssay(object = MesenchymeIntegrated.combined) <- "RNA"
nk.markers <- FindConservedMarkers(object = MesenchymeIntegrated.combined, ident.1 = 11, grouping.var = "stage", 
                                   verbose = FALSE)
head(x = nk.markers)
write.csv(nk.markers, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/Integrated/Integrated_ConservedMarkers.csv")
FeaturePlot(object = MesenchymeIntegrated.combined, features = c("Hic1","Runx2","Tfap2b","Meox2", "Shox2", "Top2a", "Dlx5", "Sp7", "Tubb3", "Wnt16", "Krt14", "Msx1", "Sox9", "Pax9", "Wwp2", "Ugdh", "Susd5", "Sox10", "Piezo1"), min.cutoff = "q9")


##### Make UMAPs with marker genes #####-------------------------
Idents(MesenchymeIntegrated.combined) <- "stage"

stage_markers <- FindAllMarkers(MesenchymeIntegrated.combined, only.pos = T, 
                                min.pct = 0.25, logfc.threshold = log2(1.5), base = 2) 
  # ^ this will identify markers which show at least 1.5 fold change between stages 

top50_stage_markers <- stage_markers %>% group_by(cluster) %>% top_n(50)
  # ^ select the top 50 genes that distinguish the stages 

MesenchymeIntegrated.combined <- RunUMAP(MesenchymeIntegrated.combined, features = top50_stage_markers$gene, 
                                         reduction.name="mumap", reduction.key="MUMAP_")
  # ^ this will make a UMAP using only the genes that are in the top 50 differential expressed genes.
  # ^ one could also select top 50 genes for each stage instead of top 50 overall and see how that affects
  # ^ the UMAP and clustering. 
  # ^ i also assign a new umap label to distinguish it from the old umap based on all variable genes.

DimPlot(MesenchymeIntegrated.combined, reduction = "mumap", group.by = "stage")


#End of integration analysis
  
  
