# Roman Thomas // 6.20.2025
# This is Homework 1, where we will be analyzing data to discuss endothelial cells 

# Change your working directory to the ROC#2, which is inside of the BTG -> Coding -> Week 1 -> Materials and Resources by clicking on "Session" (top of the screen), then go to "Set working directory", and select "Choose working directory". 
# Then choose the "desktop as your working directory 

# Activate the R packages needed to run the analysis.
# Highlight lines 8-10 by clicking on the line number 22 and dragging down to the number 24.
# To run the clode, press Ctrl+Enter (or Cmd+Enter) simultaneously to run the code.

library(Seurat)
library(dplyr)
library(patchwork)

# 1. What are the top 5 genes that define the endothelial cell population? 
# a) Create 2 different kinds of plots that would demonstrate this. 

Glom_merged <- readRDS("Glom_merged.RDS")

#double check correct load of RDS 
DimPlot(Glom_merged, reduction = "umap")

VlnPlot(Glom_merged, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln"))

#Re-name the clusters by cell-type, and add that to the metadata

new.cluster.ids <- c("Podo1", "Podo2", "Endo", "Tubule", "Mural")
names(new.cluster.ids) <- levels(Glom_merged)
Glom_merged <- RenameIdents(Glom_merged, new.cluster.ids)
DimPlot(Glom_merged, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE) + NoLegend()

#identify the top 5 genes in each cell type(cluster)
Glom_merged[["RNA"]] <- JoinLayers(Glom_merged[["RNA"]])
Glom.markers <- FindAllMarkers(Glom_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <-Glom.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5

DoHeatmap(Glom_merged, features = top5$gene) 

#The top 5 genes that define endothelial cell population are: 
#Emcn, Ctla2a, Pi16, Ramp2, and Cyp4b1

VlnPlot(Glom_merged, features = c("Emcn","Ctla2a","Pi16","Ramp2","Cyp4b1"))
#this was second plot to demonstrate 

# 2. Add a new metadata column that combines the replicate information and the cell type data. 
# use this to ensure the correct idents are present(all the cell types)
levels(Idents(Glom_merged))

Glom_merged@meta.data$cell_rep <- paste0(Idents(Glom_merged), "_", Glom_merged@meta.data$replicate)
head(Glom_merged@meta.data)

# 3. Create a table that shows the number of endothelial cells in each sample by replicate 
cell_counts_by_replicate <- table(Glom_merged@meta.data$cell_rep)
print(cell_counts_by_replicate)
# There are 615 endothelial cells for rep1 and 298 endothelial cells for rep2

# 4. Subset and re-cluster the endothelial population.
# a. How many additional clusters are there? 
# b. Are any of them contamination of other cell types? 
# c. Name the clusters based on your best guess on identity. 
# d. Create a graph that demonstrates the contaminating cell types. 

Glom_Endo <-subset(Glom_merged, idents = c("Endo"))
# use this to ensure this new object only has endothelial cells
levels(Idents(Glom_Endo))

# Re-cluster the cells in Glom_Endo, starting with data normalization
Glom_Endo <- NormalizeData(Glom_Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Glom_Endo <- FindVariableFeatures(Glom_Endo, selection.method = "vst", nfeatures = 2000)

# Now re-scale the data
# Note that the row names are now based on the row names of the new Seurat object

Endo.genes <- rownames(Glom_Endo)
Glom_tubule <- ScaleData(Glom_Endo, features = Endo.genes)

Glom_Endo <- RunPCA(Glom_Endo, features = VariableFeatures(object = Glom_Endo))

print(Glom_Endo[["pca"]], dims = 1:5, nfeatures = 5)

#Visualize the genes driving the first two PCAs
VizDimLoadings(Glom_Endo, dims = 1:2, reduction = "pca")

#Visualize the PCA
DimPlot(Glom_Endo, reduction = "pca")

#This line generates the UMAP visualization 
DimPlot(Glom_Endo, reduction = "umap")
DimPlot(Glom_merged, reduction = "umap")
# appears that there are podocyte(podo 1 and podo 2) cell types that are contaminating 

VlnPlot(Glom_Endo, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln"))
# shows that podocytes are responsible for contamination,, very prevalent in cluster 2

#Bonus: Create a UMAP that chooses colors based on the national parks color palate 
install.packages("NatParksPalettes")
library(NatParksPalettes)
library(ggplot2)
ls("package:NatParksPalettes")
class(NatParksPalettes)

print(NatParksPalettes)

umap_data <- data.frame(UMAP1 = rnorm(100), UMAP2 = rnorm(100), group = factor(rep(1:5, each = 20)))
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = group)) +
geom_point(size = 3) +  # Scatter plot with points
scale_color_manual(values = NatParksPalettes) +  # Apply NatParksPalettes color palette
labs(title = "UMAP Plot with NatParks Colors", x = "UMAP1", y = "UMAP2") + 
theme_minimal()  # Apply minimal theme

length((unique(umap_data$group)))
length(NatParksPalettes)


NatParksPalettes <- CreateSeuratObject(counts = my_list$counts, meta.data = my_list$meta_data)


DimPlot(NatParksPalettes, reduction = "umap")
