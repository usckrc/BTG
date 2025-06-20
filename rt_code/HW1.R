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
Glom_merged <- AddMetaData(object = Glom_merged, metadata = c("Podo1", "Podo2", "Endo", "Tubule", "Mural"), col.name = "cell type")
head(Glom_merged@meta.data)

# 3. Create a table that shows the number of endothelial cells in each sample by replicate 
cell_counts_by_replicate <- table(Glom_Endo@meta.data$replicate)
print(cell_counts_by_replicate)

# 4. Subset and re-cluster the endothelial population.
# a. How many additional clusters are there? 
# b. Are any of them contamination of other cell types? 
# c. Name the clusters based on your best guess on identity. 
# d. Create a graph that demonstrates the contaminating cell types. 

Glom_Endo <-subset(Glom_Endo, idents = c("Endo"))
Glom_Endo <- NormalizeData(Glom_Endo)
Glom_Endo <- FindVariableFeatures(Glom_Endo)

# Now re-scale the data
# Note that the row names are now based on the row names of the new Seurat object

Endo.genes <- rownames(Glom_Endo)
Glom_tubule <- ScaleData(Glom_Endo, features = Endo.genes)

Glom_tubule <- RunPCA(Glom_tubule, features = VariableFeatures(object = Glom_tubule))

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
names(NatParksPalettes)
