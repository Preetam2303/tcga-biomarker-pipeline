# Load libraries
library(readr)
library(dplyr)
library(biomaRt)
############
# Load the data (modify the path as needed)
file_path <- "C:/Users/KIIT/Downloads/oral cancer/Tumor.csv"
df <- read_csv(file_path)

# Rename the first column to 'GeneNames'
df <- df %>% rename(GeneNames = `...1`)
# Remove version numbers from Ensembl IDs
df$GeneNames <- sub("\\..*", "", df$GeneNames)
head(df)
# Connect to Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Ensembl IDs to gene names
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = df$GeneNames,
                      mart = mart)
##########################################################
# Merge the gene mapping with the original data
df_mapped <- merge(gene_mapping, df, by.x = "ensembl_gene_id", by.y = "GeneNames", all.y = TRUE)

# Replace missing gene names with Ensembl IDs if not found
df_mapped$hgnc_symbol[is.na(df_mapped$hgnc_symbol)] <- df_mapped$ensembl_gene_id

# Remove the 'ensembl_gene_id' column (optional)
df_mapped <- dplyr::select(df_mapped, -ensembl_gene_id)
####################
# Check dimensions
original_dim <- dim(df)
converted_dim <- dim(df_mapped)

# Print the dimensions
cat("Original data dimensions:", original_dim, "\n")
cat("Converted data dimensions:", converted_dim, "\n")
head(df_mapped)
#########################
# Remove duplicates based on 'hgnc_symbol'
df_mapped <- df_mapped %>% distinct(hgnc_symbol, .keep_all = TRUE)
# Remove rows with any missing values
df_mapped <- na.omit(df_mapped)
# Rename the 'hgnc_symbol' column to 'GeneNames'
colnames(df_mapped)[1] <- "GeneNames"
# Remove duplicate 'GeneNames' columns
df_mapped <- df_mapped[, !duplicated(names(df_mapped))]
# Copy df_mapped to df
df <- df_mapped

# Check dimensions again
converted_dim <- dim(df)
cat("Updated converted data dimensions:", converted_dim, "\n")
##################################################################
# Load libraries
library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)

# Check for missing or duplicate gene names
dim_before <- dim(df)
# Remove rows with missing or empty gene names and remove duplicates
df_filtered <- df %>%
  filter(!is.na(GeneNames) & GeneNames != "") %>%
  distinct(GeneNames, .keep_all = TRUE) %>%
  as.data.frame()  # Convert back to data frame if needed

# Check dimensions
dim_after <- dim(df_filtered)
cat("Dimensions after filtering:", dim_after, "\n")

# Check dimensions after filtering
dim_after <- dim(df_filtered)
cat("Dimensions before filtering:", dim_before, "\n")
cat("Dimensions after filtering:", dim_after, "\n")

# Set GeneNames as rownames and remove the GeneNames column
df_filtered <- as.data.frame(df_filtered)  # Convert to data frame to allow rownames
rownames(df_filtered) <- df_filtered$GeneNames
df_filtered <- df_filtered[, -1]
###########################
# Define the file path
file_path <- "C:/Users/KIIT/Downloads/oral cancer/Tumornames.csv"

# Save the filtered data to the specified file path
write.csv(df_filtered, file = file_path, row.names = TRUE)

cat("Filtered data saved to:", file_path, "\n")
# Count the number of Ensembl IDs in the row names (assuming they start with "ENS")
ensembl_count <- sum(grepl("^ENS", rownames(df_filtered)))

# Display the count of Ensembl IDs
cat("Number of Ensembl IDs in the data:", ensembl_count, "\n")

##################################
# Convert to matrix
expr_matrix <- as.matrix(df_filtered)

# Check for NA values
num_na <- sum(is.na(expr_matrix))
cat("Number of missing (NA) values in the matrix:", num_na, "\n")

# Stop execution if row names are missing
if (is.null(rownames(expr_matrix))) {
  stop("Row names (gene names) are missing from the matrix.")
}

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = expr_matrix, 
                                 project = "OralCancerRNAseq", 
                                 min.cells = 3, 
                                 min.features = 200)

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
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Visualize clusters
DimPlot(seurat_obj, reduction = "umap",pt.size = 2, label = TRUE)
saveRDS(seurat_obj, file = "C:/Users/KIIT/Downloads/oral_cancer_seurat.rds")

################################################
head(rownames(seurat_obj))
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
# Extract module scores
immune_scores <- seurat_obj@meta.data$Immune_Score1
nonimmune_scores <- seurat_obj@meta.data$NonImmune_Score1

# Classify clusters based on the higher module score
seurat_obj@meta.data$Immune_Status <- ifelse(immune_scores > nonimmune_scores, 
                                             "Immune Activation", 
                                             "Immune Evasion")

# Visualize the classification
DimPlot(seurat_obj, group.by = "Immune_Status", reduction = "umap", label = TRUE, pt.size = 2) + 
  ggtitle("Immune Status of Clusters")
###################
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
# Extract cluster information and immune scores
cluster_scores <- seurat_obj@meta.data %>%
  dplyr::select(seurat_clusters, Immune_Score1, NonImmune_Score1) %>%
  group_by(seurat_clusters) %>%
  summarise(
    Mean_Immune_Activation = mean(Immune_Score1, na.rm = TRUE),
    Mean_Immune_Evasion = mean(NonImmune_Score1, na.rm = TRUE)
  )

# Print cluster-wise scores
print(cluster_scores)

#####################################
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

#########################################
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
#############################
# Load required libraries
library(Seurat)
library(limma)

# Extract normalized expression data using 'layer' instead of 'slot'
expr_matrix <- as.matrix(GetAssayData(seurat_obj, layer = "data"))  # Log-normalized counts

# Get cluster identities
cluster_info <- Idents(seurat_obj)

# Define groups
immune_activation_clusters <- c("0", "1", "2", "3")
immune_evasion_clusters <- c("4", "5")

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
####################################
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)  # For non-overlapping gene labels

# Assume 'deg_results' contains differential expression results from limma
# Ensure it has columns: 'logFC' (log fold change) and 'adj.P.Val' (adjusted p-value)

# Define significance thresholds
logFC_threshold <- 1  # Adjust as needed
pval_threshold <- 0.05

# Assign DEG results to limma_results
limma_results <- deg_results

# Ensure gene names are present
limma_results$gene <- rownames(limma_results)  

# Create a new column for gene classification
limma_results <- limma_results %>%
  mutate(GeneCategory = case_when(
    adj.P.Val < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
    adj.P.Val < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 10 most significant DEGs (smallest adjusted p-value)
top_degs <- limma_results %>%
  filter(GeneCategory != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  head(10)  # Select top 10 most significant DEGs

# Define colors for the plot
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")

# Plot volcano plot with top 10 gene labels
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.7, size = 2) +  # Increase point size for better visibility
  geom_text_repel(data = top_degs, aes(label = gene), size = 3, box.padding = 0.3, max.overlaps = 10) +  # Label top 10 genes
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value",
       color = "DEG Category") +  # Custom legend title
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "gray", size = 1) +
  theme(legend.position = "right")
##########################################
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define significance thresholds
logFC_threshold <- 1  
pval_threshold <- 0.05

# Ensure gene names are present
limma_results$gene <- rownames(limma_results)  

# Classify genes
limma_results <- limma_results %>%
  mutate(GeneCategory = case_when(
    adj.P.Val < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
    adj.P.Val < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 10 most significant DEGs
top_degs <- limma_results %>%
  filter(GeneCategory != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  head(10)

# Define colors
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")

# Create volcano plot
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.6, size = 1) +  # Smaller points for better clarity
  geom_text_repel(data = top_degs, aes(label = gene), size = 3, box.padding = 0.5, max.overlaps = 10) +  
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "logFC",
       y = "-log10(adj.P.Val)",
       color = "DEG Category") +  
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(legend.position = "right")
####################################
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define significance thresholds
logFC_threshold <- 0.05 
pval_threshold <- 0.05

# Ensure gene names are present
limma_results$gene <- rownames(limma_results)  

# Avoid infinite -log10(adj.P.Val) by setting a small lower limit
limma_results$adj.P.Val <- pmax(limma_results$adj.P.Val, 1e-100)

# Classify genes
limma_results <- limma_results %>%
  mutate(GeneCategory = case_when(
    adj.P.Val < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
    adj.P.Val < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 5 most significant DEGs in each category for labeling
top_up <- limma_results %>% filter(GeneCategory == "Upregulated") %>% arrange(adj.P.Val) %>% head(10)
top_down <- limma_results %>% filter(GeneCategory == "Downregulated") %>% arrange(adj.P.Val) %>% head(10)
top_degs <- bind_rows(top_up, top_down)

# Define colors
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")

# Create volcano plot
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.6, size = 1) +  # Smaller points for better clarity
  geom_text_repel(data = top_degs, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3) +  
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "logFC",
       y = "-log10(adj.P.Val)",
       color = "DEG Category") +  
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(legend.position = "right") +
  ylim(0, 50)  # Adjust Y-axis range to avoid excessive stretching
#################################
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.6, size = 1) +  # Smaller points for better clarity
  geom_text_repel(data = top_degs, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3) +  
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "logFC",
       y = "-log10(adj.P.Val)",
       color = "DEG Category") +  
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(legend.position = "right") +
  ylim(0, 50) +  # Adjust Y-axis
  xlim(-0.95,0.75)    # Adjust X-axis
############################
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define significance thresholds
logFC_threshold <- 0.05 
pval_threshold <- 0.05

# Ensure gene names are present
limma_results$gene <- rownames(limma_results)  

# Avoid infinite -log10(adj.P.Val) by setting a small lower limit
limma_results$adj.P.Val <- pmax(limma_results$adj.P.Val, 1e-100)

# Classify genes
limma_results <- limma_results %>%
  mutate(GeneCategory = case_when(
    adj.P.Val < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
    adj.P.Val < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 10 most significant DEGs in each category for labeling
top_up <- limma_results %>% filter(GeneCategory == "Upregulated") %>% arrange(adj.P.Val) %>% head(10)
top_down <- limma_results %>% filter(GeneCategory == "Downregulated") %>% arrange(adj.P.Val) %>% head(10)
top_degs <- bind_rows(top_up, top_down)

# Extract top 10 most significant miRNAs (genes starting with "MIR")
top_mirna <- limma_results %>% filter(grepl("^MIR", gene)) %>% arrange(adj.P.Val) %>% head(10)

# Define colors
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")

# Create volcano plot
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.6, size = 1) +  # Smaller points for better clarity
  # Label top DEGs (Upregulated & Downregulated)
  geom_text_repel(data = top_degs, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3) +  
  # Label miRNAs in GREEN
  geom_text_repel(data = top_mirna, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3, color = "green") +  
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "logFC",
       y = "-log10(adj.P.Val)",
       color = "DEG Category") +  
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(legend.position = "right") +
  ylim(0, 50) +  # Adjust Y-axis range to avoid excessive stretching
  xlim(-0.95,0.75)    # Adjust X-axis
####################################
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define significance thresholds
logFC_threshold <- 0.05 
pval_threshold <- 0.05

# Ensure gene names are present
limma_results$gene <- rownames(limma_results)  

# Avoid infinite -log10(adj.P.Val) by setting a small lower limit
limma_results$adj.P.Val <- pmax(limma_results$adj.P.Val, 1e-100)

# Classify genes
limma_results <- limma_results %>%
  mutate(GeneCategory = case_when(
    adj.P.Val < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
    adj.P.Val < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 10 most significant DEGs in each category for labeling
top_up <- limma_results %>% filter(GeneCategory == "Upregulated") %>% arrange(adj.P.Val) %>% head(10)
top_down <- limma_results %>% filter(GeneCategory == "Downregulated") %>% arrange(adj.P.Val) %>% head(10)
top_degs <- bind_rows(top_up, top_down)

# Extract top 10 most significant miRNAs (genes starting with "MIR")
top_mirna <- limma_results %>% filter(grepl("^MIR", gene)) %>% arrange(adj.P.Val) %>% head(10)

# Define colors
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black", "miRNA" = "green")

# Assign miRNAs a separate category for coloring points
limma_results <- limma_results %>%
  mutate(GeneCategory = ifelse(gene %in% top_mirna$gene, "miRNA", GeneCategory))

# Create volcano plot
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.6, size = 1) +  # Smaller points for better clarity
  # Label top DEGs (Upregulated & Downregulated)
  geom_text_repel(data = top_degs, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3) +  
  # Label miRNAs in BLACK
  geom_text_repel(data = top_mirna, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3, color = "black") +  
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "logFC",
       y = "-log10(adj.P.Val)",
       color = "DEG Category") +  
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(legend.position = "right") +
  ylim(0, 50) +  # Adjust Y-axis range to avoid excessive stretching
  xlim(-0.95,0.75)    # Adjust X-axis
############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define significance thresholds
logFC_threshold <- 0.05 
pval_threshold <- 0.05

# Ensure gene names are present
limma_results$gene <- rownames(limma_results)  

# Avoid infinite -log10(adj.P.Val) by setting a small lower limit
limma_results$adj.P.Val <- pmax(limma_results$adj.P.Val, 1e-100)

# Explicitly create GeneCategory column
limma_results <- limma_results %>%
  mutate(GeneCategory = case_when(
    grepl("^MIR", gene) ~ "miRNA",  # Assign all miRNAs to the "miRNA" category
    adj.P.Val < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
    adj.P.Val < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Check if GeneCategory was created successfully
print(table(limma_results$GeneCategory))  # Debugging step

# Select top 10 most significant DEGs in each category for labeling
top_up <- limma_results %>% filter(GeneCategory == "Upregulated") %>% arrange(adj.P.Val) %>% head(10)
top_down <- limma_results %>% filter(GeneCategory == "Downregulated") %>% arrange(adj.P.Val) %>% head(10)
top_degs <- bind_rows(top_up, top_down)

# Extract top 10 most significant miRNAs
top_mirna <- limma_results %>% filter(GeneCategory == "miRNA") %>% arrange(adj.P.Val) %>% head(10)

# Define colors
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black", "miRNA" = "green")

# Create volcano plot
ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = GeneCategory)) +
  geom_point(alpha = 0.6, size = ifelse(limma_results$GeneCategory == "miRNA", 4, 1)) +  # Larger points for miRNAs
  # Label top DEGs (Upregulated & Downregulated)
  geom_text_repel(data = top_degs, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3) +  
  # Label miRNAs in BLACK
  geom_text_repel(data = top_mirna, aes(label = gene), 
                  size = 3, box.padding = 0.3, max.overlaps = Inf, segment.size = 0.3, color = "black") +  
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "logFC",
       y = "-log10(adj.P.Val)",
       color = "DEG Category") +  
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(legend.position = "right") +
  ylim(0, 50) +  # Adjust Y-axis range to avoid excessive stretching
  xlim(-0.95,0.75)    # Adjust X-axis
############################################
head(limma_results)
str(limma_results)
table(limma_results$GeneCategory)  # Check if it exists
print(top_mirna_deg)
##########
# Extract expression matrix for DEGs only
deg_expression <- expr_matrix[rownames(degs_filtered), , drop = FALSE]  

# Scale expression values across samples (Z-score transformation)
deg_expression_scaled <- t(scale(t(deg_expression)))  

# Check dimensions to ensure correct subset
dim(deg_expression_scaled)
###########################
library(pheatmap)

pheatmap(
  deg_expression_scaled,  
  scale = "none",  # Already scaled
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  color = colorRampPalette(c("blue", "white", "red"))(50),  
  main = "Heatmap of Differentially Expressed Genes (DEGs)"
)
#######################
sample_groups <- data.frame(Group = factor(group_labels))
rownames(sample_groups) <- colnames(deg_expression_scaled)  # Match column names

pheatmap(
  deg_expression_scaled,  
  scale = "none",  
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  annotation_col = sample_groups,  # Add group labels  
  color = colorRampPalette(c("blue", "white", "red"))(50),  
  main = "Heatmap of Differentially Expressed Genes (DEGs)"
)
#################
pheatmap(
  deg_expression_scaled,  
  scale = "none",  
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  annotation_col = sample_groups,  
  color = colorRampPalette(c("blue", "white", "red"))(50),  
  main = "Heatmap of Differentially Expressed Genes (DEGs)",  
  show_colnames = FALSE  # Hides sample names  
)
##################################
head(seurat_obj@meta.data)
# Install ComplexHeatmap from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# Install circlize from CRAN
install.packages("circlize")

###################
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Get immune scores from Seurat object
immune_scores <- seurat_obj@meta.data$Immune_Score1
nonim_scores <- seurat_obj@meta.data$NonImmune_Score1

# Calculate the difference between immune and non-immune scores
immune_status <- immune_scores - nonim_scores

# Order samples based on immune status
sample_order <- order(immune_status, decreasing = TRUE)

# Get the genes from degs_filtered
deg_genes <- rownames(degs_filtered)

# Extract expression matrix for these DEGs
expr_matrix <- GetAssayData(seurat_obj, slot = "data")[deg_genes, ]

# Scale the expression matrix
expr_scaled <- t(scale(t(expr_matrix)))

# Create annotation for immune status
ha_column = HeatmapAnnotation(
  Immune_Status = immune_status,
  col = list(Immune_Status = colorRamp2(
    c(min(immune_status), 0, max(immune_status)),
    c("blue", "white", "red")
  )),
  show_legend = TRUE,
  annotation_name_side = "left"
)

# Create color function for heatmap
col_fun = colorRamp2(
  c(-2, 0, 2), 
  c("blue", "white", "red")
)

# Create the heatmap
heatmap <- Heatmap(
  expr_scaled,
  name = "Expression",
  col = col_fun,
  cluster_columns = FALSE,  # Don't cluster columns as we want to preserve the immune score order
  cluster_rows = TRUE,      # Cluster rows to show gene patterns
  show_column_names = FALSE,
  show_row_names = FALSE,   # Hide row names due to large number of DEGs
  column_title = "Samples ordered by Immune Status",
  row_title = "Differentially Expressed Genes",
  top_annotation = ha_column,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# Draw the heatmap
draw(heatmap, padding = unit(c(2, 2, 2, 2), "mm"))

# Save the plot if needed
pdf("immune_status_DEGs_heatmap.pdf", width = 12, height = 10)
draw(heatmap, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# If you want to show only top DEGs (e.g., top 50 by absolute log fold change)
top_degs <- degs_filtered %>%
  arrange(desc(abs(logFC))) %>%
  head(50) %>%
  rownames()

# Create heatmap with only top DEGs
expr_matrix_top <- GetAssayData(seurat_obj, slot = "data")[top_degs, ]
expr_scaled_top <- t(scale(t(expr_matrix_top)))

# Create heatmap for top DEGs
heatmap_top <- Heatmap(
  expr_scaled_top,
  name = "Expression",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,   # Show row names for top DEGs
  column_title = "Samples ordered by Immune Status",
  row_title = "Top 50 Differentially Expressed Genes",
  top_annotation = ha_column,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# Draw the heatmap with top DEGs
draw(heatmap_top, padding = unit(c(2, 2, 2, 2), "mm"))

# Save the top DEGs plot
pdf("immune_status_top_DEGs_heatmap.pdf", width = 12, height = 10)
draw(heatmap_top, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
#################################################
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Get immune scores from Seurat object
immune_scores <- seurat_obj@meta.data$Immune_Score1
nonim_scores <- seurat_obj@meta.data$NonImmune_Score1

# Calculate the difference between immune and non-immune scores
immune_status <- immune_scores - nonim_scores

# Order samples based on immune status (ascending order)
sample_order <- order(immune_status)  # Removed decreasing=TRUE to get ascending order

# Get top 50 DEGs by absolute log fold change
top_degs <- degs_filtered %>%
  arrange(desc(abs(logFC))) %>%
  head(50) %>%
  rownames()

# Extract expression matrix for top DEGs
expr_matrix_top <- GetAssayData(seurat_obj, slot = "data")[top_degs, ]
expr_scaled_top <- t(scale(t(expr_matrix_top)))

# Create annotation for immune status
ha_column = HeatmapAnnotation(
  Immune_Status = immune_status,
  col = list(Immune_Status = colorRamp2(
    c(min(immune_status), 0, max(immune_status)),
    c("blue", "white", "red")
  )),
  show_legend = TRUE,
  annotation_name_side = "left"
)

# Create color function for heatmap
col_fun = colorRamp2(
  c(-2, 0, 2), 
  c("blue", "white", "red")
)

# Create heatmap for top DEGs
heatmap_top <- Heatmap(
  expr_scaled_top[, sample_order],  # Use ordered samples
  name = "Expression",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_title = "Samples ordered by Immune Status (Low to High)",
  row_title = "Top 50 Differentially Expressed Genes",
  top_annotation = ha_column,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# Draw the heatmap
draw(heatmap_top, padding = unit(c(2, 2, 2, 2), "mm"))

# Save the plot
pdf("immune_status_top_DEGs_heatmap_ascending.pdf", width = 12, height = 10)
draw(heatmap_top, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
