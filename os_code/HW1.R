########################################################
#Merged Data creation
########################################################

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# We will be uploading the data slightly differently this time, by using the read.table function 
# Convert the .txt file into a data file

Glom_rep1.data <- read.table(file = "GSM3022239_dge_glom_rep1.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric")))

# Checkpoint: 22827 genes and 4623 cells 

#########################
#########################
#header = TRUE, Indicates that the first row of the file contains column names.  
#row.names = 1, Uses the first column in the file as the row names (e.g., gene IDs).
#########################
#########################


#Glom_rep1.data is a data frame that can be utilized in several different ways
#We will next create a seurat object out of the data frame
Glom_rep1 <- CreateSeuratObject(counts = Glom_rep1.data, project = "Glom_merged", min.cells = 3, min.features = 200)

Glom_rep1

#########################
#########################
#•	counts = Glom_rep1.data
#Supplies the gene expression matrix (likely raw UMI counts) as input. Rows = genes, columns = cells.
#•	project = "Glom_merged"
#Labels this dataset/project as "Glom_merged"—useful for downstream tracking.
#•	min.cells = 3
#Filters out genes that are expressed in fewer than 3 cells.
#•	min.features = 200
#Filters out cells that have fewer than 200 detected genes. These are likely poor-quality cells or empty droplets.
#########################
#########################

#An object of class Seurat 
#17511 features across 4556 samples within 1 assay 
#Active assay: RNA (17511 features, 0 variable features)

#Since we are planning to merge this dataset with another, we will need to check to see
#if its identity is coded in the metadata

head(Glom_rep1@meta.data)

#As you can see, the column orig.ident contains "Glom_merged", which is is the project name
#Once we add another dataset, that will be overwritten, and can't be used to pull out data
#We will add a column to the metadata that denotes the replicate

Glom_rep1 <- AddMetaData(object = Glom_rep1, metadata = "rep1", col.name = "replicate")
head(Glom_rep1@meta.data)

#########################
#########################
#•	AddMetaData(...)
#This Seurat function adds new metadata to each cell in your Seurat object.
#•	object = Glom_rep1
#The Seurat object you’re adding metadata to.
#•	metadata = "rep1"
#assigning the string "rep1" as metadata. Since it’s a single value, it will be added to every cell in the object.
#•	col.name = "replicate"
#The name of the new metadata column will be "replicate".
#########################
#########################

#Next, make the second second replicate a data frame
Glom_rep2.data <- read.table(file = "GSM3022240_dge_glom_rep2.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric")))

#Glom_rep1.data is a data frame that can be utilized in several different ways
#We will next create a seurat object out of the dataframe
Glom_rep2 <- CreateSeuratObject(counts = Glom_rep2.data, project = "Glom_merged", min.cells = 3, min.features = 200)

Glom_rep2

#An object of class Seurat 
#18309 features across 4366 samples within 1 assay 
#Active assay: RNA (18309 features, 0 variable features)

#Add a column to the metadata indicate the replicate
Glom_rep2 <- AddMetaData(object = Glom_rep2, metadata = "rep2", col.name = "replicate")
head(Glom_rep2@meta.data)

# Merge the two Seurat objects. We will use the "merge" function, but there are several different ways to merge Seurat objects
# Note: the "merge" function does not work when datasets are very different from one another - i.e., stim cell vs. unstim cell
# It does seem to work when comparing disease vs. no disease
########################################################
########################################################
#?????????????How can I merge more than two data sets?????????????????????
########################################################
########################################################
Glom_merged <- merge(Glom_rep1, y = Glom_rep2, add.cell.ids = c("Glom_rep1", "Glom_rep2"), project = "Glom_merged")
Glom_merged

########################################################
########################################################
#add.cell.ids = c("Glom_rep1", "Glom_rep2")
#•	This adds a prefix to each cell’s name from each dataset, like:
#•	Glom_rep1_AAAC... for cells from Glom_rep1
#•	Glom_rep2_TTCC... for cells from Glom_rep2
########################################################
########################################################

#An object of class Seurat 
#19183 features across 8922 samples within 1 assay 
#Active assay: RNA (19183 features, 0 variable features)
#Note that the sample number is the sum of the two reps: 4556 + 4366 = 8922

# notice the cell names now have an added identifier
head(colnames(Glom_merged))

# We can also double check the number of cells in each group using the "table" function, pulling out "replicate"
table(Glom_merged$replicate)

# Now, begin the QC that we did last month, starting with filtering out samples with a high % of mito genes
# This adds a column in the metadata called "percent.mt"
Glom_merged[["percent.mt"]] <- PercentageFeatureSet(Glom_merged, pattern = "^mt-")
head(Glom_merged@meta.data)

########################################################
########################################################
#PercentageFeatureSet(Glom_merged, pattern = "^mt-")
#•	This function calculates the percentage of mitochondrial gene expression for each cell.
#•	pattern = "^mt-" tells Seurat to find all genes starting with “mt-”, which is the common #prefix for mitochondrial genes (e.g., mt-Nd1, mt-Cytb in mouse).
#•	It sums up all mitochondrial gene counts for each cell and divides that by the total gene #counts for the same cell, returning a percentage.

#Glom_merged[["percent.mt"]] <- ...
#•	This stores the output (a numeric vector of mitochondrial percentages per cell) into the #metadata slot of your Seurat object.
#•	percent.mt becomes a new column in Glom_merged@meta.data, so you can use it later for:
#•	Filtering out low-quality cells
#•	QC plots (e.g., violin plots of percent.mt vs nFeature_RNA)
########################################################
########################################################


# As with last month visualize the dataset for quality control using 3 different metrics: 
#   nFeature_RNA = Number of differing genes detected in each cell
#   nCount_RNA = Number of mRNA molecules detected in each cell
#   percent.mt = Percent of counts that come from mitochondrial genes

VlnPlot(Glom_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "replicate")
VlnPlot(Glom_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
plots <- lapply(features, function(feat) {
  data <- FetchData(Glom_merged, vars = c(feat, "replicate"))
  
  # Assign custom alpha values
  data$alpha <- ifelse(data$replicate == "rep2", 0.3, 1.0)
  
  ggplot(data, aes(x = "", y = .data[[feat]], color = replicate, alpha = alpha)) +
    geom_violin(fill = "gray90", width = 1.2, color = NA) +  # gray violin, no border
    geom_jitter(width = 0.2, size = 0.5) +
    scale_color_manual(values = c("rep1" = "black", "rep2" = "red")) +
    guides(alpha = "none") +  # hide alpha from legend
    theme_minimal() +
    labs(y = feat, x = "") +
    theme(legend.position = "right")
})
wrap_plots(plots, ncol = 3)


########################################################
########################################################
#🔢 features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
#•	Defines the list of features to plot (as strings)
#•	These are metadata metrics per cell:
#•	nFeature_RNA: number of genes detected
#•	nCount_RNA: total RNA counts
#•	percent.mt: percent of mitochondrial gene expression
#🔁 plots <- lapply(features, function(feat) { ... })
#•	Loops through each feature in the features list
#•	For each, generates a custom violin plot and stores it in a list
#📤 data <- FetchData(Glom_merged, vars = c(feat, "replicate"))
#•	Extracts a data frame from the Seurat object containing:
#•	The current feature’s values
#•	The replicate column (used to color points)
#🎨 data$alpha <- ifelse(data$replicate == "rep2", 0.3, 1.0)
#•	Adds a new alpha column to control point transparency:
#•	Cells from "rep2" get opacity of 0.3
#•	All others (e.g. "rep1") get full opacity (1.0)
#📊 ggplot(data, aes(...)) + ...
#•	Builds the violin plot using ggplot2
#aes(x = "", y = .data[[feat]], color = replicate, alpha = alpha)
#•	All cells go into a single category on the x-axis
#•	y-axis shows the current feature’s values
#•	Color is mapped to replicate (for black/red coloring)
#•	Alpha (opacity) is based on the new alpha column
#geom_violin(fill = "gray90", width = 1.2, color = NA)
#•	Draws a light gray violin (distribution shape)
#•	No border line (color = NA)
#geom_jitter(width = 0.2, size = 0.5)
#•	Plots each cell as a point with jitter to avoid overlap
#•	Uses the color and alpha mappings defined above
#scale_color_manual(values = c("rep1" = "black", "rep2" = "red"))
#•	Manually sets the point colors:
#•	"rep1" → black
#•	"rep2" → red
#guides(alpha = "none")
#•	Hides the alpha (transparency) legend
#theme_minimal()
#•	Applies a clean, minimal plot theme
#labs(y = feat, x = "")
#•	Sets the y-axis label to the feature name
#•	Removes x-axis label
#theme(legend.position = "right")
#•	Moves the color legend to the right side of the plot
#🧱 wrap_plots(plots, ncol = 3)
#•	Combines all the individual plots into a single 3-column layout using patchwork
########################################################
########################################################


# We will interpret this information to filter out doublets and damaged cells
#   High nFeature_RNA or nCount_RNA could equate to doublets or mixed cell debris. We will filter out cells with more than 2000 genes and 4000 mRNA molecules
#   High percent.mt could equate to damaged cells. We will filter out cells with more than 20% mitochondrial content.

Glom_merged <- subset(Glom_merged, subset = nFeature_RNA < 2000 & nCount_RNA < 4000 & percent.mt < 20)

# No we will replot the graphs to visualize the plots again with the newly filtered data

VlnPlot(Glom_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Glom_merged
#An object of class Seurat 
#19183 features across 8361 samples within 1 assay 
#Active assay: RNA (19183 features, 0 variable features)

#This filtered out 561 cells and 0 genes

# Data Normalization and Identification of Genes Driving Variability (cluster formation)
# Normalize the matrix and create a plot of the genes that are the most variable in the dataset

Glom_merged <- NormalizeData(Glom_merged, normalization.method = "LogNormalize", scale.factor = 10000)

Glom_merged <- FindVariableFeatures(Glom_merged, selection.method = "vst", nfeatures = 2000)

########################################################
########################################################
# 🔹 Glom_merged <- NormalizeData(Glom_merged, normalization.method = "LogNormalize", scale.factor = 10000)
# •	Goal: Make gene expression values comparable across cells
# •	What it does for each cell:
# •	Add up all gene counts
# •	Divide each gene by the total → this standardizes across cells
# •	Multiply by 10,000 → puts values on a similar scale
# •	Take the log of each value → reduces effect of very large numbers
# 🔹 Glom_merged <- FindVariableFeatures(...)
# •	Identifies the top 2,000 most variable genes
# •	Uses the “vst” method (variance-stabilizing transformation)
# •	These genes will be used in downstream steps (like PCA)
########################################################
########################################################


# This line Identifies the 10 most highly variable genes and creates a new object in the enviornment called "top10"
top10 <- head(VariableFeatures(Glom_merged), 10)

# These lines creates a plot of the variable features within the dataset and labels the top 10 genes to see which are driving the clustering
plot1 <- VariableFeaturePlot(Glom_merged)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Creation and Visualization of Principal Components
# Calculate the principal components and then print the 5 genes that drive the first 5 principal componsents.
# Can you get a sense for how the cells might be clustering based off of these genes?
# Note: the "ScaleData" and "RunPCA" steps can take very long with large datasets

all.genes <- rownames(Glom_merged)
Glom_merged <- ScaleData(Glom_merged, features = all.genes)

Glom_merged <- RunPCA(Glom_merged, features = VariableFeatures(object = Glom_merged))

print(Glom_merged[["pca"]], dims = 1:5, nfeatures = 5)

########################################################
########################################################
# 🔹 all.genes <- rownames(Glom_merged)
# •	Retrieves all gene names from the Seurat object Glom_merged
# •	rownames(Glom_merged) gives the row names of the RNA assay — which are gene symbols
# •	Stores them in a variable called all.genes for later use (e.g., in scaling)
# 🔹 Glom_merged <- ScaleData(Glom_merged, features = all.genes)
# •	Performs z-score scaling on the gene expression values for the listed features (all.genes)
# •	For each gene:
# •	Subtracts the mean expression across all cells
# •	Divides by the standard deviation
# •	This step is essential before PCA because:
# •	PCA assumes all features (genes) are on a comparable scale
# •	Without this, highly expressed genes would dominate the PCs just due to scale
# •	This does not change raw data, just creates a scaled version used internally
# 🔹 Glom_merged <- RunPCA(Glom_merged, features = VariableFeatures(object = Glom_merged))
# •	Runs Principal Component Analysis (PCA) on the scaled data
# •	Uses only the variable genes (not all genes) — these were identified earlier using FindVariableFeatures()
# •	PCA finds new axes (principal components) that explain the most variance in the data
# •	Each PC is a linear combination of genes, and each cell gets a score along each PC
# •	The result is saved to Glom_merged[["pca"]], which stores:
#   •	PC scores per cell (for clustering, UMAP, etc.)
# •	Loadings per gene (to see which genes drive each PC)
# 🔹 print(Glom_merged[["pca"]], dims = 1:5, nfeatures = 5)
# •	Prints out PCA gene loadings (i.e., top genes contributing to each PC)
# •	For PC1 to PC5, it shows the top 5 genes that have the highest weights (“loadings”)
# •	These are the genes most responsible for the separation of cells along those PCs
# •	Helpful for:
# •	Understanding what biological signals each PC might represent
# •	Checking which cell types or processes may be driving variation
########################################################
########################################################


#Now we will visualize the genes driving the PCA a few ways


#here we can visualize the PCA scores in a table.
pca_scores <- Embeddings(Glom_merged, reduction = "pca")
View(pca_scores)
pca_scores

#This plot shows the variability and helps to decide how many PCA dimensions to use when measuring the clusters. 
#When the standard deviation in the PC dimensions becomes small it no longer has a big effect. In this case we will use 10 PC's.
ElbowPlot(Glom_merged)

########################################################
########################################################
########################################################
########################################################
#Extract standard deviations of PCs
stdev <- Glom_merged[["pca"]]@stdev
#Compute % variance explained
pct_var <- (stdev^2) / sum(stdev^2) * 100
#Make a data frame
pc_df <- data.frame(
  PC = factor(paste0("PC", 1:length(pct_var)), 
  levels = paste0("PC", 1:length(pct_var))),
  Variance = pct_var,
  Label = sprintf("%.1f%%", pct_var)
  )
#Plot as a bar plot
ggplot(pc_df[1:10, ], aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Label), vjust = -0.3, size = 3.5) +  # label above bar
  labs(title = "Scree Plot (Bar)", y = "% Variance Explained", x = "Principal Components") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.1))

########################################################
########################################################
# 🔹 stdev <- Glom_merged[["pca"]]@stdev
# •	Accesses the PCA results stored in the Seurat object Glom_merged
# •	["pca"]: selects the PCA reduction
# •	@stdev: extracts the standard deviation for each principal component (PC)
# •	Each standard deviation reflects how much variation that PC explains
# 🔹 pct_var <- (stdev^2) / sum(stdev^2) * 100
# •	Calculates the percentage of total variance explained by each PC
# •	stdev^2: squares each standard deviation to get the variance per PC
# •	sum(stdev^2): computes the total variance across all PCs
# •	(stdev^2) / sum(...) * 100: gives the percent contribution of each PC to the total variance
# 🔹 pc_f <- data.frame(...)
# •	Creates a data frame called pc_f to store information for plotting or analysis
# 🔸 PC = factor(paste0("PC", 1:length(pct_var)))
# •	Creates PC labels like "PC1", "PC2", …, "PC30" depending on the number of PCs
# •	paste0() joins “PC” with the index numbers
# •	factor() converts the labels to a categorical variable, which helps preserve order in plots
# 🔸 Variance = pct_var
# •	Adds a column containing the percent variance explained by each PC
# •	Matches each "PC" label with its corresponding variance
# 🔸 Label = sprintf("%.1f%%", pct_var)
# •	Formats the variance values into percent labels like "12.3%"
# •	sprintf("%.1f%%", ...):
#   •	Rounds to 1 decimal place
# •	Adds a % symbol
# •	Useful for labeling bars in a scree plot
########################################################
########################################################

########################################################
########################################################
########################################################
########################################################

#This line visualizes the genes driving the first 2 PCA dimensions
VizDimLoadings(Glom_merged, dims = 1:2, reduction = "pca")

#This line visualizes the scRNAseq dataset as a PCA (hopefully you can appreciate the poor clustering)
DimPlot(Glom_merged, reduction = "pca", group.by = "replicate")
DimPlot(Glom_merged, reduction = "pca")

########################################################
########################################################
# •	DimPlot(): a Seurat function to visualize dimensionality reductions
# •	Glom_merged: your Seurat object containing PCA results
# •	reduction = "pca": tells Seurat to use the PCA coordinates for plotting
########################################################
########################################################

#This line creates a heatmap of the top 500 genes that drive the first 2 PCA dimensions. 
DimHeatmap(Glom_merged, dims = 1:2, cells = 500, balanced = TRUE)


# Calculation and Creation of a UMAP visualization
# Determine the shape of the clusters and the number of populations in our scRNASeq visualization.

#These lines control the resolution of the clusters and the number of separate populations.
# "FindNeighbors" changes the separation between clusters
# "FindClusters" changes the number of differntial populations
# DISCUSSION ABOUT 10 for line 166
Glom_merged <- FindNeighbors(Glom_merged, dims = 1:10)
Glom_merged <- FindClusters(Glom_merged, resolution = 0.15)
Glom_merged <- RunUMAP(Glom_merged, dims = 1:10)
# This line generates the UMAP visualization. How do you think the clusters look?
DimPlot(Glom_merged, reduction = "umap")

########################################################
########################################################
# 🔹 Glom_merged <- FindNeighbors(Glom_merged, dims = 1:10)
# •	Purpose: Calculates a cell-to-cell similarity graph based on PCA
# •	Function: FindNeighbors() from Seurat
# •	dims = 1:10: Uses PCs 1 through 10 for distance calculation
# •	Builds a k-nearest neighbors (KNN) graph of cells in PCA space
# •	✅ This graph is used as the foundation for clustering (next step)
# 🔹 Glom_merged <- FindClusters(Glom_merged, resolution = 0.1)
# •	Purpose: Assigns cells to clusters based on the neighbor graph
# •	Function: FindClusters() from Seurat
# •	Uses Louvain algorithm or Leiden algorithm for community detection
# •	resolution = 0.1:
# •	Controls how many clusters you get
# •	Higher value = more clusters, lower value = fewer clusters
# •	✅ This assigns each cell a cluster identity (stored in seurat_clusters)
# 🔹 Glom_merged <- RunUMAP(Glom_merged, dims = 1:10)
# •	Purpose: Performs UMAP (Uniform Manifold Approximation and Projection) for 2D visualization
# •	Function: RunUMAP() from Seurat
# •	dims = 1:10: Uses the first 10 PCs as input
# •	Reduces high-dimensional data into a 2D space that preserves structure and clusters
# •	✅ Allows you to visualize your cells in a UMAP plot
# 🔹 DimPlot(Glom_merged, reduction = "umap")
# •	Purpose: Plots the UMAP embedding
# •	Function: DimPlot() from Seurat
# •	reduction = "umap": Specifies that UMAP coordinates should be used
# •	Each point = one cell
# •	Colors = different clusters (from FindClusters)
# •	✅ Lets you see how cells are grouped visually and interpret biological structure
########################################################
########################################################

# The first sanity/QC check for merged datasets is to see if cells are clustering based on replicate
# Since we have replicate info in the metadata, we can check this using the "group.by" function
DimPlot(Glom_merged, reduction = "umap", group.by = "replicate")

# The replicates are integrated nicely. It does not appear there is any batch effect
# If there is a batch effect, consider using the "FindIntegrationAnchors" function, which we will not address in today's meeting
# https://satijalab.org/seurat/archive/v3.0/immune_alignment.html
# You can also decrease the number of variable features in the "FindVariableFeatures" function. Right now it is set to 2000, decreasing it
# will help deal with batch effect, but will decrease the resolution of clustering

# We will next go through the cell identification features we covered in last month's meeting

# First generation Heatmaps of differential genes between clusters to identify cell types
# This can take a VERY long time with large datasets

Glom_merged[["RNA"]] <- JoinLayers(Glom_merged[["RNA"]])
Glom.markers <- FindAllMarkers(Glom_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Glom.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- Glom.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(Glom_merged, features = top5$gene) 

########################################################
########################################################
# 🔹 Glom_merged[["RNA"]] <- JoinLayers(Glom_merged[["RNA"]])
# •	Purpose: Combines RNA assay layers (e.g. “counts”, “data”) into one coherent layer
# •	Function: JoinLayers() is used in Seurat v5 to merge multi-layered assay data
# •	Makes sure everything (like raw counts and normalized data) is in sync
# •	✅ Required in Seurat v5 before running FindAllMarkers()
# 🔹Glom.markers <- FindAllMarkers(Glom_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# •	Purpose: Finds marker genes for each cluster
# •	Function: FindAllMarkers() compares each cluster vs. all others
# •	Parameters:
#   •	only.pos = TRUE: only returns genes upregulated in each cluster
# •	min.pct = 0.25: includes genes expressed in at least 25% of cells in a group
# •	logfc.threshold = 0.25: filters genes with log2 fold-change > 0.25
# •	✅ Returns a table with marker genes, p-values, avg log2FC, etc.
# 🔹 Glom.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# •	Purpose: Selects the top 5 marker genes per cluster
# •	Function:
#   •	group_by(cluster): groups the markers by cluster ID
# •	top_n(n = 5, wt = avg_log2FC): selects the top 5 genes per cluster with the highest log2 fold-change
# •	✅ Allows you to highlight the most defining genes for each cluster
# 🔹 top5 <- Glom.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
# •	Same as above, but saves the result into a new variable top5
# •	top5$gene will contain the names of the most important genes per cluster
########################################################
########################################################


# Nphs2 = Podocyte; Pecam1 = Endothelial; Slc12a3 = Tubule; Tagln = Mural
VlnPlot(Glom_merged, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln"))

# It looks like there are two different podocyte populations
# We can check out there differences between clusters 0 and 1 using the "FindMarkers" function
cluster0.markers <- FindMarkers(Glom_merged, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
head(cluster0.markers, n = 25)

#Using this information, we can visualize the differences using the "Violin Plot" function
VlnPlot(Glom_merged, features = c("Ctgf"))

# We won't spend a lot of time trying to identify different podocyte populations
# But you could look at populations more distinct, such as endothelial (population 2) and epithelial (population 3) cells
cluster2.markers <- FindMarkers(Glom_merged, ident.1 = 2, ident.2 = 3, min.pct = 0.25)
head(cluster2.markers, n = 10)
VlnPlot(Glom_merged, features = c("Atp1b1"))

# Re-name the clusters by cell-type, and add that to the metadata

new.cluster.ids <- c("Podo1", "Podo2", "Endo", "Tubule", "Mural", "Immune")
names(new.cluster.ids) <- levels(Glom_merged)
Glom_merged <- RenameIdents(Glom_merged, new.cluster.ids)
DimPlot(Glom_merged, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE) + NoLegend()

# Note that Podo1 and Podo2 are right next to each other, yet the Endo population is split apart
# I personally think this means the two podocyte populations are not significant
# The other Endo "cluster" may represent doublets, but we will not spend much time addressing this
# If you want, you can spend time with the "FeaturePlot" function to visualize what genes are expressed in the other Endo "Cluster", 
# then adjust your "FindNeighbors" and "FindClusters" values a bit.

FeaturePlot(Glom_merged, features = c("Nphs2") )






########################################################
########################################################
#QUESTION 1:What are the top 5 genes that define the endothelial cell population? Create 2 different kinds of plots that would demonstrate this.
########################################################
########################################################

#STEP1
Glom_merged[["RNA"]] <- JoinLayers(Glom_merged[["RNA"]])
Glom.markers <- FindAllMarkers(Glom_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Glom.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- Glom.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(Glom_merged, features = top5$gene) 

#STEP2
Glom_Endo <-subset(Glom_merged, idents = c("Endo"))
Glom_Endo

#STEP3
top5_Endo <- Glom.markers %>% filter(cluster == 2) %>% top_n(5, avg_log2FC)
length(top5_Endo$gene)
print(top5_Endo$gene)

#STEP4
DoHeatmap(Glom_Endo, features = top5$gene) 
DoHeatmap(Glom_Endo, features = top5_Endo$gene) 

#Srgn Ctla2a Pi16 Ramp2 Cyp4b1 

#STEP5
top5_Endo_genes <- unique(top5_Endo$gene)[1:5]
VlnPlot(Glom_merged, features = c("Srgn", "Ctla2a", "Pi16", "Ramp2", "Cyp4b1"), pt.size = 0.1, ncol = 2)
VlnPlot(Glom_merged, features = c("Cdh5", "Pecam1", "Pi16", "Plat", "Cyp4b1",  "Ehd3", "Kdr", "Nr2f2", "Emcn"), pt.size = 0.1, ncol = 3)
VlnPlot(Glom_merged, features = c("Plat", "Emcn", "Tsapn7", "Kdr", "Bmx", "Mapt", "Smad6", "Ehd3", "Lpl", "Flt1", "Fbln2", "Mgp", "Trpv4"), pt.size = 0.1, ncol = 4)
 

########################################################
########################################################
#QUESTION 2:Add a new metadata column that combines the replicate information and the cell type data.
########################################################
########################################################

#STEP1 Create a metadata column for cell type from the current identities
Glom_merged$celltype <- Idents(Glom_merged)

#STEP2 Combine replicate and cell type into a new metadata column
Glom_merged$rep_celltype <- paste(Glom_merged$replicate, Glom_merged$celltype, sep = "_")
head(Glom_merged@meta.data)


########################################################
########################################################
#QUESTION 3:Create a table that shows the number of endothelial cells in each sample by replicate
########################################################
########################################################

table(Glom_merged$replicate)
table(Glom_merged$replicate, Glom_merged$celltype)

sum(Idents(Glom_merged) == "Endo")


########################################################
########################################################
#QUESTION 4:Subset and re-cluster the endothelial population.
# a.	How many additional clusters are there? 2
# b.	Are any of them contamination of other cell types?
# c.	Name the clusters based on your best guess on identity.
# d.	Create a graph that demonstrates the contaminating cell types.
########################################################
########################################################

Glom_Endo <-subset(Glom_merged, idents = c("Endo"))
Glom_Endo

Glom_Endo <- NormalizeData(Glom_Endo, normalization.method = "LogNormalize", scale.factor = 10000)

########################################################
########################################################
# 🔹 Glom_Endo <- NormalizeData(Glom_Endo, normalization.method = "LogNormalize", scale.factor = 10000)
# •	Purpose: Normalizes the gene expression values for each cell.
# •	Function: NormalizeData() scales each cell’s total gene expression to the scale.factor (10,000), then log-transforms it.
# •	Why: Accounts for differences in sequencing depth across cells.
# •	Detailed Explanation:Normalization adjusts for the fact that some cells may have more RNA captured or sequenced than others. NormalizeData() performs a per-cell correction by scaling each cell’s total counts to a common value (e.g., 10,000), then applying log-transformation. This ensures that differences in gene expression across cells reflect biological differences, not technical ones. It makes gene expression values comparable across cells.
########################################################
########################################################

Glom_Endo <- FindVariableFeatures(Glom_Endo, selection.method = "vst", nfeatures = 2000)

########################################################
########################################################
# 🔹 Glom_Endo <- FindVariableFeatures(Glom_Endo, selection.method = "vst", nfeatures = 2000)
# •	Purpose: Identifies the top 2,000 highly variable genes across endothelial cells.
# •	Function:
# •	selection.method = "vst" uses variance-stabilizing transformation to detect variable genes.
# •	nfeatures = 2000 selects the top 2000 most variable genes.
# •	Why: These genes are more informative for PCA and clustering than uniformly expressed genes.
########################################################
########################################################

endo.genes <- rownames(Glom_Endo)

########################################################
########################################################
# 🔹 endo.genes <- rownames(Glom_Endo)
# •	Purpose: Saves all gene names (rows) in Glom_Endo to a variable.
# •	Why: You’ll use this list to scale all genes (not just variable ones).
########################################################
########################################################

Glom_Endo <- ScaleData(Glom_Endo, features = endo.genes)

########################################################
########################################################
# 🔹 Glom_Endo <- ScaleData(Glom_Endo, features = endo.genes)
# •	Purpose: Centers and scales the expression data for each gene (z-score).
# •	Function: ScaleData() subtracts the mean and divides by the standard deviation across all cells.
# •	Why: Necessary for PCA to work properly, since genes must be on a comparable scale.
# •	Detailed Explanation:While NormalizeData() corrects for total RNA content between cells, ScaleData() standardizes gene expression across all cells for each gene. It centers each gene’s expression around 0 and scales it to unit variance. This ensures that all genes contribute equally to downstream dimensionality reduction (like PCA). Without this, genes with higher variance would dominate the analysis, even if that variance isn’t biologically meaningful. Together, NormalizeData() and ScaleData() make your data comparable both across cells and across genes.
########################################################
########################################################

Glom_Endo <- RunPCA(Glom_Endo, features = VariableFeatures(object = Glom_Endo))
print(Glom_Endo[["pca"]], dims = 1:5, nfeatures = 5)

########################################################
########################################################
# 🔹 Glom_Endo <- RunPCA(Glom_Endo, features = VariableFeatures(object = Glom_Endo))
#   • Purpose: Performs principal component analysis on the scaled data using the 2,000 variable genes.
#   • Function: RunPCA() reduces the dimensionality of the data to capture the major axes of variation.
# 🔹 print(Glom_Endo[["pca"]], dims = 1:5, nfeatures = 5)
#   • Purpose: Prints the top 5 genes that contribute to each of the first 5 principal components.
#   • Function: Helps you see which genes are driving variation in each PC.
########################################################
########################################################

VizDimLoadings(Glom_Endo, dims = 1:2, reduction = "pca")

########################################################
########################################################
# 🔹 VizDimLoadings(Glom_Endo, dims = 1:2, reduction = "pca")
#   • Purpose: Visualizes gene loadings (contributions) to PCs 1 and 2.
#   • Result: Bar plots of the genes contributing most to each PC.
########################################################
########################################################

DimPlot(Glom_Endo, reduction = "pca")

########################################################
########################################################
# 🔹 DimPlot(Glom_Endo, reduction = "pca")
#   • Purpose: Plots cells in PC1 vs. PC2 space.
#   • Why: Allows visual inspection of how cells separate based on gene expression.
########################################################
########################################################

DimHeatmap(Glom_Endo, dims = 1:2, cells = 500, balanced = TRUE)

########################################################
########################################################
# 🔹 DimHeatmap(Glom_Endo, dims = 1:2, cells = 500, balanced = TRUE)
#   • Purpose: Heatmap of top genes contributing to PCs 1 and 2 across a subset of 500 cells.
#   • balanced = TRUE: Ensures equal number of up- and down-regulated genes per PC are shown.
#   • Use: Helps visually inspect biological meaning of PCs.
########################################################
########################################################

ElbowPlot(Glom_Endo)

########################################################
########################################################
# 🔹 ElbowPlot(Glom_Endo)
#   • Purpose: Shows the variance explained by each PC (scree plot).
#   • How to use: Helps decide how many PCs to use for clustering/UMAP by finding the “elbow” point where additional PCs add little value.
########################################################
########################################################

Glom_Endo <- FindNeighbors(Glom_Endo, dims = 1:7)

########################################################
########################################################
# 🔹 Glom_Endo <- FindNeighbors(Glom_Endo, dims = 1:7)
#   • Purpose: Builds a nearest-neighbor graph based on the first 7 PCs.
#   • Function: Used for clustering and UMAP.
#   • Why 7 PCs? Based on elbow plot judgment.
########################################################
########################################################

Glom_Endo <- FindClusters(Glom_Endo, resolution = 0.1)

########################################################
########################################################
# 🔹 Glom_Endo <- FindClusters(Glom_Endo, resolution = 0.1)
#   • Purpose: Performs Louvain clustering on the neighbor graph.
#   • resolution = 0.1: Controls how many clusters you get — lower = fewer clusters.
########################################################
########################################################

Glom_Endo <- RunUMAP(Glom_Endo, dims = 1:7)

########################################################
########################################################
# 🔹 Glom_Endo <- RunUMAP(Glom_Endo, dims = 1:7)
#   • Purpose: Runs UMAP dimensionality reduction to visualize clusters in 2D.
#   • Input: Uses first 7 PCs to determine cell similarity.
########################################################
########################################################

DimPlot(Glom_Endo, reduction = "umap")
DimPlot(Glom_Endo, reduction = "umap", group.by = "replicate")

length(unique(Idents(Glom_Endo))) #how many clusters
table(Idents(Glom_Endo))#how many cells per cluster

# Let's do a sanity check just to make sure we aren't doing the same cluster all over again
# Do a violin plot of the same genes that defined the original clusteres
VlnPlot(Glom_Endo, features = c("Nrp1", "Cdh5", "Eln", "Plat", "Emcn", "Bmx", "Trpv4", "Mapt", "Fbln2", "Edn1", "Cryab", "Kdr"))

# We will double check our work by visualizing some "traditional" cell-type markers with violin plots by running line 155
# Nphs2 = Podocyte; Pecam1 = Endothelial; Slc12a3 = Tubule; Tagln = Mural
VlnPlot(Glom_Endo, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln"))

#Next we can visualize the expression of each gene on the UMAP diagram by running the following lines individually.         
FeaturePlot(Glom_Endo, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln","Ptprc" ),label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE)


#these are the glomerular endo markers
FeaturePlot(Glom_Endo, features = c("Plat"),label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE)
FeaturePlot(Glom_Endo, features = c("Emcn"),label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE)


#subsetting again from the endo cluster 0 subset to achieve pure results
Glom_Endo_clean <- subset(Glom_Endo, idents = "0")
Glom_Endo_clean <- NormalizeData(Glom_Endo_clean)
Glom_Endo_clean <- FindVariableFeatures(Glom_Endo_clean)
Glom_Endo_clean <- ScaleData(Glom_Endo_clean)
Glom_Endo_clean <- RunPCA(Glom_Endo_clean)
ElbowPlot(Glom_Endo_clean)
Glom_Endo_clean <- FindNeighbors(Glom_Endo_clean, dims = 1:10)
Glom_Endo_clean <- FindClusters(Glom_Endo_clean, resolution = 0.1)  # Slightly higher resolution to catch fine states
Glom_Endo_clean <- RunUMAP(Glom_Endo_clean, dims = 1:10)
DimPlot(Glom_Endo_clean, label = TRUE)

Glom_Endo_clean[["RNA"]] <- JoinLayers(Glom_Endo_clean[["RNA"]])
Glom_Endo_clean.markers <- FindAllMarkers(Glom_Endo_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Glom_Endo_clean.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_Glom_Endo_clean <- Glom_Endo_clean.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(Glom_Endo_clean, features = top5_Glom_Endo_clean$gene) 

VlnPlot(Glom_Endo_clean, features = c("Plat", "Emcn", "Ehd3", "Fbln2", "Fbln5", "Cldn5", "Efnb2", "Mgp", "Edn1", "Tspan8"))
FeaturePlot(Glom_Endo_clean, features = c("Plat", "Emcn", "Ehd3", "Fbln2", "Fbln5", "Cldn5", "Mgp"),label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE)

new.ids <- c("G", "A or VA")
names(new.ids) <- levels(Glom_Endo_clean)
Glom_Endo_clean <- RenameIdents(Glom_Endo_clean, new.ids)
DimPlot(Glom_Endo_clean, label = TRUE)
