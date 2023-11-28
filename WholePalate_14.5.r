---
  title: "E14.5 whole palate"
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
E14.5whole.data <- Read10X(data.dir = "C:/Users/siddh/OneDrive - University of Southern California/Palate/filtered_feature_bc_matrix/E14.5/filtered_feature_bc_matrix/")

#create seurat object
E14.5whole <- CreateSeuratObject(counts = E14.5whole.data, project = "E14.5_WT")
```

```{r}
E14.5whole[["percent.mt"]] <- PercentageFeatureSet(E14.5whole, pattern = "^mt-")

VlnPlot(E14.5whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


```{r}

FeatureScatter(E14.5whole, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(E14.5whole, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


dim(E14.5whole)#genes=32285 cells=5023
sum(E14.5whole$percent.mt < 15)

```

```{r}
E14.5whole <- subset(E14.5whole, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 15)
E14.5whole <- NormalizeData(E14.5whole, normalization.method = "LogNormalize", scale.factor = 10000)
E14.5whole <- FindVariableFeatures(E14.5whole, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(E14.5whole)
E14.5whole <- ScaleData(E14.5whole, features = all.genes)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(E14.5whole), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(E14.5whole)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2



```

```{r}
#load ccgenes 
"/Single cell/Soft palate E14.5 by Summer/E14.5/Mouse cell cycle genes/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds"
mouse_cell_cycle_genes <- readRDS("C:/Users/siddh/OneDrive - University of Southern California/Palate/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes[[1]]
g2m.genes <- mouse_cell_cycle_genes[[2]]
```

```{r}
E14.5whole <- CellCycleScoring(E14.5whole, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E14.5whole <- RunPCA(E14.5whole, features = c(s.genes, g2m.genes))
p1 <- DimPlot(E14.5whole)
E14.5whole <- ScaleData(E14.5whole, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E14.5whole))
E14.5whole <- RunPCA(E14.5whole, features = VariableFeatures(E14.5whole), nfeatures.print = 10)
E14.5whole <- RunPCA(E14.5whole, features = c(s.genes, g2m.genes))
p2 <- DimPlot(E14.5whole)

p1+p2
```

```{r}
E14.5whole <- RunPCA(E14.5whole, features = VariableFeatures(object = E14.5whole))
ElbowPlot(E14.5whole)
```

```{r}
E14.5whole <- FindNeighbors(E14.5whole, dims = 1:30)
E14.5whole <- FindClusters(E14.5whole, resolution = 0.5)
```

```{r}
E14.5whole <- RunUMAP(E14.5whole, dims = 1:0)
```

```{r}
DimPlot(E14.5whole, reduction = "umap", label = T, label.size = 4, group.by = "RNA_snn_res.0.5")
```
```{r}
#
saveRDS(E14.5whole, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E14.5/E14.5whole_afterCC.rds")
```
E14.5whole <- readRDS(file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E14.5/E14.5whole_afterCC.rds")

```{r}
DimPlot(E14.5whole, reduction = "umap", group.by = "Phase")
```


```{r}
FeaturePlot(E14.5whole, features = c("Tfap2b", "Bglap", "Runx2", "Sp7"), order = T)
```

```{r}
FeaturePlot(E14.5whole, features = c("Tfap2b", "Bglap", "Runx2", "Sp7"), order = T)
```
FeaturePlot(E14.5whole, features = c("Aldh1a2","Hic1","Tbx15","Fgf18","Tfap2b","Tbx22", "Meox1", "Dlx5", "Sp7","Top2a"), order = T)
FeaturePlot(E14.5whole, features = c("Krt14","Tubb3","Hba-x","Lyz2","Sox10","Cdh5", "Dlx5", "Top2a", "Twist1", "Piezo2"), order = T)
# Perimysial markers- "Aldh1a2","Hic1","Tbx15"
# CNC marker= "Tfap2b"
# Mitotic cells= "Top2a"
# CNC derived mesenchymal cells= "Meox1", "Dlx5"


```{r}
FeaturePlot(E14.5whole, features = c("Runx2", "Sp7"), blend = T, blend.threshold = 0, cols= c("blue", "red"))
```

```{r}
E14.5whole.markers <- FindAllMarkers(E14.5whole, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_E14.5whole <- E14.5whole.markers %>% group_by(cluster) %>% top_n(n = 10)
top100_E14.5whole <- E14.5whole.markers %>% group_by(cluster) %>% top_n(n = 100)
DoHeatmap(E14.5whole, features = top10_E14.5whole$gene) 

```
```{r}
#
write.csv(top10_E14.5whole, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E14.5/top10_E14.5.whole.csv")
```

```
Mesenchyme_E14.5 <- subset(E13.5whole, idents = c("1","10","4","0","2"))
DimPlot(Mesenchyme_E14.5, reduction = "umap", label = TRUE, pt.size = 0.5) 
saveRDS(Mesenchyme_E14.5, file = "C:/Users/siddh/OneDrive - University of Southern California/Palate/Analysis/E14.5/Dim20/Mesenchyme_E14.5.rds")
```
