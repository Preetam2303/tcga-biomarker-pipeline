# Load necessary library
library(Seurat)

# Load the Seurat object
seurat_obj <- readRDS(file = r"(C:\Users\KIIT\Downloads\oral cancer\immune activation vs immune evasion\oral_cancer_seurat.rds)")

# Verify the object
print(seurat_obj)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

# Identify highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)

# Scale the data
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:12)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:16,min.dist = 0.5,spread = 0.9)

# Visualize clusters
DimPlot(seurat_obj, reduction = "umap",pt.size = 2, label = TRUE)
########################

########################
immune_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2",
                  "ERAP1", "ERAP2", "TAP1", "TAP2", "NLRC5", "IRF1", "NFKB1", "RELA", "CXCL9",
                  "CXCL10", "CCL2", "CCL5", "MB21D1", "TMEM173", "IFI16", "TBK1", "IRF3", "HLA-E")

non_immune_genes <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CD244", "CD160", "IRF4", "BATF", 
                      "NFATC1", "FOXP3", "IL2RA", "TNFRSF18", "IKZF4", "CD80", "CD86", "CD274", 
                      "PDCD1LG2", "LGALS9")
# Check immune genes
present_immune_genes <- intersect(immune_genes, rownames(seurat_obj))
cat("Immune genes present in the dataset:\n")
print(present_immune_genes)

# Check non-immune genes
present_non_immune_genes <- intersect(non_immune_genes, rownames(seurat_obj))
cat("\nNon-immune genes present in the dataset:\n")
print(present_non_immune_genes)
####################################################################
# Calculate immune module score
seurat_obj <- AddModuleScore(seurat_obj, features = list(immune_genes), name = "Immune_Score")

# Calculate non-immune module score
seurat_obj <- AddModuleScore(seurat_obj, features = list(non_immune_genes), name = "NonImmune_Score")

# Visualize the module scores using a violin plot
VlnPlot(seurat_obj, features = c("Immune_Score1", "NonImmune_Score1"), pt.size = 0.1, ncol = 2)
#########################
# Add cluster information to the metadata
cluster_metadata <- data.frame(Cluster = Idents(seurat_obj))

# Combine cluster metadata with Immune and Non-Immune scores
seurat_obj$Cluster <- cluster_metadata$Cluster

# Calculate average scores for Immune and Non-Immune genes across clusters
immune_scores_avg <- aggregate(Immune_Score1 ~ Cluster, data = seurat_obj[[]], mean)
nonimmune_scores_avg <- aggregate(NonImmune_Score1 ~ Cluster, data = seurat_obj[[]], mean)

# Print the average scores
cat("Average Immune Scores by Cluster:\n")
print(immune_scores_avg)

cat("\nAverage Non-Immune Scores by Cluster:\n")
print(nonimmune_scores_avg)

##############################3
# Extract module scores
immune_scores <- seurat_obj@meta.data$Immune_Score1
nonimmune_scores <- seurat_obj@meta.data$NonImmune_Score1

# Classify clusters based on the higher module score
seurat_obj@meta.data$Immune_Status <- ifelse(immune_scores > nonimmune_scores, 
                                             "Immune Activation", 
                                             "Immune Evasion")
library(ggplot2)
library(patchwork)
# Visualize the classification
DimPlot(seurat_obj, group.by = "Immune_Status", reduction = "umap", label = TRUE, pt.size = 2) + 
  ggtitle("Immune Status of Clusters")
##############################
# Load required libraries
library(Seurat)
library(ggplot2)

# Define cluster colors
cluster_colors <- c("0" = "red", "1" = "blue", "2" = "green", 
                    "3" = "purple", "4" = "orange", "5" = "yellow")

# Create a data frame for custom legend labels
legend_labels <- data.frame(
  cluster = factor(c("0", "1", "2", "3", "4", "5"), levels = c("0", "1", "2", "3", "4", "5")),
  immune_status = c("Immune Activation", "Immune Activation", "Immune Activation", "Immune Activation",
                    "Immune Evasion", "Immune Evasion")
)

# Create the UMAP plot
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", 
                     cols = cluster_colors, pt.size = 2) + 
  labs(title = "Immune Classification of Clusters") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12, face = "bold")) +
  guides(color = guide_legend(title = "Immune Status", 
                              override.aes = list(size = 3)))  # Decrease color symbol size

# Modify the legend to display clusters under the correct immune categories
umap_plot + 
  scale_color_manual(values = cluster_colors, 
                     breaks = legend_labels$cluster,
                     labels = paste0(legend_labels$immune_status, "\n", legend_labels$cluster)) +
  theme(legend.text = element_text(size = 10))
########################################
library(ggplot2)
library(Seurat)

# Define immune classifications
immune_activation_clusters <- c("0", "1", "2", "3")
immune_evasion_clusters <- c("4", "5")

# Assign classification to Seurat metadata
seurat_obj@meta.data$Immune_Classification <- ifelse(
  seurat_obj$seurat_clusters %in% immune_activation_clusters, "Immune Activation",
  "Immune Evasion"
)

# Define custom colors
immune_colors <- c("Immune Activation" = "red", "Immune Evasion" = "blue")

# Generate UMAP plot with immune classification
DimPlot(seurat_obj, group.by = "Immune_Classification", reduction = "umap", label = FALSE, pt.size = 2) + 
  scale_color_manual(values = immune_colors) +
  ggtitle("Immune Classification of Clusters") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  labs(color = "Immune Status\n(Red - Immune Activation: Clusters 0-3\nBlue - Immune Evasion: Clusters 4,5)")
#################################
# Load required libraries
library(Seurat)
library(limma)

# Extract normalized expression data using 'layer' instead of 'slot'
expr_matrix <- as.matrix(GetAssayData(seurat_obj, layer = "data"))  # Log-normalized counts

# Get cluster identities
cluster_info <- Idents(seurat_obj)

# Define groups
immune_activation_clusters <- c("0", "1", "2")
immune_evasion_clusters <- c("3", "4")

# Create sample group labels
group_labels <- ifelse(cluster_info %in% immune_activation_clusters, "Activation", "Evasion")

# Create design matrix
design <- model.matrix(~ 0 + factor(group_labels))
colnames(design) <- c("Activation", "Evasion")

# Fit the linear model using limma
fit <- lmFit(expr_matrix, design)

# Define contrast (Activation vs Evasion)
contrast_matrix <- makeContrasts(Activation - Evasion, levels = design)

# Apply contrast to the model
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract differentially expressed genes (DEGs)
deg_results <- topTable(fit2, adjust.method = "fdr", number = Inf)

# Filter DEGs based on logFC and FDR
logFC_threshold <- 0.75  # Change as needed
FDR_threshold <- 0.05
degs_filtered <- deg_results[abs(deg_results$logFC) > logFC_threshold & deg_results$adj.P.Val < FDR_threshold, ]

# View the top DEGs
head(degs_filtered)
dim(degs_filtered)  # Get the number of rows and columns

# Save results to a CSV file
write.csv(degs_filtered, "DEGs_Activation_vs_Evasion.csv")
############################################################
# Load required libraries
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)

# Load the DEGs data
degs_data <- read.csv("C:/Users/KIIT/Downloads/oral cancer/immune activation vs immune evasion/DEGs_Activation_vs_Evasion.csv", row.names = 1)

# Extract the cluster identity of each cell
cluster_info <- Idents(seurat_obj)

# Extract expression data for DEGs from the Seurat object
degs_genes <- rownames(degs_data)  # Assuming row names contain gene names
expr_matrix <- GetAssayData(seurat_obj, slot = "data")[degs_genes, ]

# Scale the expression data (for better visualization)
expr_matrix_scaled <- t(scale(t(as.matrix(expr_matrix))))  # Scale across genes

# Create annotation for clusters
annotation_col <- data.frame(Cluster = factor(cluster_info))
rownames(annotation_col) <- colnames(expr_matrix_scaled)

# Plot heatmap
pheatmap(expr_matrix_scaled, 
         annotation_col = annotation_col, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row",
         color = colorRampPalette(c("purple", "black", "yellow"))(50))
##############################
# Order columns by cluster
ordered_cells <- order(annotation_col$Cluster)
expr_matrix_sorted <- expr_matrix_scaled[, ordered_cells]
annotation_col_sorted <- annotation_col[ordered_cells, , drop = FALSE]

# Plot heatmap
pheatmap(expr_matrix_sorted, 
         annotation_col = annotation_col_sorted, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Prevent clustering of columns
         scale = "row",
         show_rownames = FALSE,  # Remove row names
         color = colorRampPalette(c("purple", "black", "yellow"))(50))
##################################
# Order columns by cluster
ordered_cells <- order(annotation_col$Cluster)
expr_matrix_sorted <- expr_matrix_scaled[, ordered_cells]
annotation_col_sorted <- annotation_col[ordered_cells, , drop = FALSE]

# Get cluster labels
column_labels <- annotation_col_sorted$Cluster

# Plot heatmap
pheatmap(expr_matrix_sorted, 
         annotation_col = annotation_col_sorted, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Maintain cluster order
         scale = "row",
         show_rownames = TRUE,  # Keep row names
         show_colnames = FALSE, # Remove column names
         labels_col = column_labels,  # Mark clusters with numbers above heatmap
         color = colorRampPalette(c("purple", "black", "yellow"))(50))

#######################################
# Select top differentially expressed genes
top_degs <- head(rownames(seurat_obj), 50)  # Modify based on your DEGs list

# Generate heatmap
DoHeatmap(seurat_obj, features = top_degs, group.colors = NULL) + 
  scale_fill_gradientn(colors = c("purple", "black", "yellow"))
##########################################
# Identify markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
library(dplyr)
library(ggrepel)
# Select the top 20 genes per cluster
top20_genes <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>% 
  pull(gene)

# Generate heatmap using these top genes
DoHeatmap(seurat_obj, features = top20_genes, group.by = "seurat_clusters") + 
  scale_fill_gradientn(colors = c("purple", "black", "yellow"))
####################################
# Filter genes that start with "MIR"
mir_genes <- top20_genes[grep("^MIR", top20_genes)]

# Print the MIR genes (if any)
print(mir_genes)
###########################Don't
# Count the number of ENSEMBL gene IDs in the markers dataset
ensembl_count <- sum(grepl("^ENSG", markers$gene))

# Print the count
print(ensembl_count)
####################

#######
library(ggplot2)
library(patchwork)  # For combining plots

# Sample Data (replace this with your actual data)
set.seed(123)
umap_data_filtered <- data.frame(
  umap_1 = rnorm(300),
  umap_2 = rnorm(300),
  IA_score = rnorm(300, mean = 5, sd = 1),  # Immune Activation Score
  IE_score = rnorm(300, mean = 1, sd = 1)   # Immune Evasion Score
)

# Plot 1: Immune Activation Score
p1 <- ggplot(umap_data_filtered, aes(x = umap_1, y = umap_2, color = IA_score)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "blue", limits = c(-2, 8))+
  labs(title = "Immune Activation Score", color = "IA_Score") +
  theme_minimal()

# Plot 2: Immune Evasion Score
p2 <- ggplot(umap_data_filtered, aes(x = umap_1, y = umap_2, color = IE_score)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "blue", limits = c(-2, 8))+
  labs(title = "Immune Evasion Score", color = "IE_Score") +
  theme_minimal()

# Combine the plots side by side
p1 + p2
############################
# Add immune scores to metadata
seurat_obj@meta.data$Mean_Immune_Activation <- seurat_obj@meta.data$Immune_Score1
seurat_obj@meta.data$Mean_Immune_Evasion <- seurat_obj@meta.data$NonImmune_Score1
library(ggplot2)

# UMAP plot colored by Immune Activation score
FeaturePlot(seurat_obj, features = "Immune_Score1", reduction = "umap", cols = c("lightgrey", "blue")) + 
  ggtitle("Immune Activation Score across Clusters")

# UMAP plot colored by Immune Evasion score
FeaturePlot(seurat_obj, features = "NonImmune_Score1", reduction = "umap", cols = c("lightgrey", "red")) + 
  ggtitle("Immune Evasion Score across Clusters")
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 2) + 
  ggtitle("UMAP Clustering with Labels")
library(patchwork)

p1 <- FeaturePlot(seurat_obj, features = "Immune_Score1", reduction = "umap", cols = c("lightgrey", "blue"))
p2 <- FeaturePlot(seurat_obj, features = "NonImmune_Score1", reduction = "umap", cols = c("lightgrey", "red"))

p1 + p2
p1 <- FeaturePlot(seurat_obj, features = "Immune_Score1", reduction = "umap", cols = c("lightgrey", "blue")) + 
  ggtitle("Immune Activation")

p2 <- FeaturePlot(seurat_obj, features = "NonImmune_Score1", reduction = "umap", cols = c("lightgrey", "blue")) + 
  ggtitle("Immune Evasion")

# Combine the plots
library(patchwork)
p1 + p2
#################################
library(Seurat)
library(ggplot2)
library(patchwork)

# Feature plot for Immune Activation
p1 <- FeaturePlot(seurat_obj, features = "Immune_Score1", reduction = "umap", cols = c("lightgrey", "blue")) + 
  ggtitle("Immune Activation")

# Feature plot for Immune Evasion
p2 <- FeaturePlot(seurat_obj, features = "NonImmune_Score1", reduction = "umap", cols = c("lightgrey", "blue")) + 
  ggtitle("Immune Evasion")

# Add cluster labels using DimPlot (UMAP with clusters)
cluster_labels <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)  

# Arrange in a better layout
(p1 + cluster_labels) / (p2 + cluster_labels) + plot_layout(heights = c(1, 1))
(p1 | cluster_labels) / (p2 | cluster_labels)
