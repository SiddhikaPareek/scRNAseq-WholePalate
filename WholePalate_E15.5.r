---
  title: "E15.5 whole palate"
author: "EJ"
date: "09/19/2020"
output: html_document
---
  ```{r}
library(Seurat)
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(patchwork)
```

```{r}
#load data library no1
E15.5whole.data <- Read10X(data.dir = "C:/Users/siddh/OneDrive - University of Southern California/Palate/filtered_feature_bc_matrix/E15.5/filtered_feature_bc_matrix/")

#create seurat object
E15.5whole <- CreateSeuratObject(counts = E15.5whole.data, project = "E15.5_WT")
```

```{r}
E15.5whole[["percent.mt"]] <- PercentageFeatureSet(E15.5whole, pattern = "^mt-")

VlnPlot(E15.5whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



```{r}


FeatureScatter(E15.5whole, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(E15.5whole, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dim(E15.5whole)#genes=32285 cells=4789
sum(E15.5whole$percent.mt < 15)


```{r}
E15.5whole <- subset(E15.5whole, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 15)
E15.5whole <- NormalizeData(E15.5whole, normalization.method = "LogNormalize", scale.factor = 10000)
E15.5whole <- FindVariableFeatures(E15.5whole, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(E15.5whole)
E15.5whole <- ScaleData(E15.5whole, features = all.genes)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(E15.5whole), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(E15.5whole)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2



```

```{r}
#load ccgenes 
"/Single cell/Soft palate E15.5 by Summer/E15.5/Mouse cell cycle genes/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds"
mouse_cell_cycle_genes <- readRDS("C:/Users/siddh/OneDrive - University of Southern California/Palate/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes[[1]]
g2m.genes <- mouse_cell_cycle_genes[[2]]
```

```{r}
E15.5whole <- CellCycleScoring(E15.5whole, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E15.5whole <- RunPCA(E15.5whole, features = c(s.genes, g2m.genes))
p1 <- DimPlot(E15.5whole)
E15.5whole <- ScaleData(E15.5whole, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E15.5whole))
E15.5whole <- RunPCA(E15.5whole, features = VariableFeatures(E15.5whole), nfeatures.print = 10)
E15.5whole <- RunPCA(E15.5whole, features = c(s.genes, g2m.genes))
p2 <- DimPlot(E15.5whole)

p1+p2
```

```{r}
E15.5whole <- RunPCA(E15.5whole, features = VariableFeatures(object = E15.5whole))
ElbowPlot(E15.5whole)
```

```{r}
E15.5whole <- FindNeighbors(E15.5whole, dims = 1:20)
E15.5whole <- FindClusters(E15.5whole, resolution = 0.5)
```

```{r}
E15.5whole <- RunUMAP(E15.5whole, dims = 1:20)
```

```{r}
DimPlot(E15.5whole, reduction = "umap", label = T, label.size = 4, group.by = "RNA_snn_res.0.5")
```
```{r}
#
saveRDS(E15.5whole, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E15.5/E15.5whole_afterCC.rds")
```
E15.5whole <- readRDS(file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E15.5/E15.5whole_afterCC.rds")

```{r}
DimPlot(E15.5whole, reduction = "umap", group.by = "Phase")
```

```{r}
FeaturePlot(E15.5whole, features = c("Tfap2b", "Bglap", "Runx2", "Sp7"), order = T)
```

```{r}
FeaturePlot(E15.5whole, features = c("Tfap2b", "Bglap", "Runx2", "Sp7"), order = T)
```
FeaturePlot(E14.5whole, features = c("Tfap2b", "Bglap", "Runx2", "Sp7"), order = T)
```
FeaturePlot(E15.5whole, features = c("Aldh1a2","Hic1","Tbx15","Fgf18","Tfap2b","Tbx22", "Meox1", "Dlx5", "Sp7","Top2a"), order = T)
FeaturePlot(E15.5whole, features = c("Krt14","Tubb3","Hba-x","Lyz2","Sox10","Cdh5", "Dlx5", "Top2a", "Twist1", "Piezo2"), order = T)
# Perimysial markers- "Aldh1a2","Hic1","Tbx15"
# CNC marker= "Tfap2b"
# Mitotic cells= "Top2a"
# CNC derived mesenchymal cells= "Meox1", "Dlx5"





```{r}
FeaturePlot(E15.5whole, features = c("Runx2", "Sp7"), blend = T, blend.threshold = 0, cols= c("blue", "red"))
```

```{r}
E15.5whole.markers <- FindAllMarkers(E15.5whole, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_E15.5whole <- E15.5whole.markers %>% group_by(cluster) %>% top_n(n = 10)
top100_E15.5whole <- E15.5whole.markers %>% group_by(cluster) %>% top_n(n = 100)
DoHeatmap(E15.5whole, features = top10_E15.5whole$gene) 

```
```{r}
#
write.csv(top10_E15.5whole, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E15.5/top10_E15.5.whole.csv")
```
```
Mesenchyme_E15.5 <- subset(E13.5whole, idents = c("2","0","1","6","7","3"))
DimPlot(Mesenchyme_E15.5, reduction = "umap", label = TRUE, pt.size = 0.5) 
saveRDS(Mesenchyme_E15.5, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E15.5/Dim20/Mesenchyme_E15.5.rds")
```