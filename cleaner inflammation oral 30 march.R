# Load necessary library
library(Seurat)
library(dplyr)
library(ggplot2)

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
# Plot the elbow plot
ElbowPlot(seurat_obj, ndims = 50)

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:12)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:16,min.dist = 0.5,spread = 0.9)

# Visualize clusters
DimPlot(seurat_obj, reduction = "umap",pt.size = 2, label = TRUE)
########################
inflammatory_genes <- c("ACSL4", "AKT", "AP1", "BCL2", "C1R", "C1S", "CCFB", "CCL", "CCL2", "CCL3",
                        "CCL4", "CCL5", "CCL7", "CCL11", "CCL27", "CD1C", "CD3", "CD4", "CD8B", "CD10", 
                        "CD11B", "CD11C", "CD14", "CD15", "CD16", "CD19", "CD20", "CD24", "CD25", 
                        "CD27", "CD45RA", "CD45RO", "CD80", "CD81", "CD86", "CD138", "CD141", "CD163", 
                        "CD200", "CD206", "CF", "CFB", "CJUN", "COX2", "CTACK", "CTLA4", "CXCL", 
                        "CXCL1", "CXCL2", "CXCL3", "CXCL6", "CXCL8", "CXCL9", "CXCL10", "CXCL12", 
                        "CXCL14", "CXCL16", "CXCR1", "CXCR2", "CXCR7", "EGFR", "FGF2", "FN1", "FOXP3", 
                        "FNI", "GARP", "GBP1", "GBPI", "HGF", "HLA-DRA", "HRAS", "ICAM1", "IFI44L", 
                        "IFIT3", "IFNG", "IL1", "IL1B", "IL1RN", "IL2", "IL2RA", "IL4", "IL5", "IL6", 
                        "IL7", "IL8", "IL9", "IL10", "IL11", "IL12", "IL12B", "IL13", "IL16", "IL17", 
                        "IL17A", "IL18", "ISG15", "KITLG", "KRTI5", "LAP3", "LIF", "LTA", "MIF", 
                        "MMP1", "MMP3", "MMP7", "MMP10", "MTARC2", "NGF", "NFKB", "NMI", "OASL", 
                        "PDGFA", "PDGFB", "PD-1", "PDL-1", "PGE2", "PLAU", "SCFD1", "SCGB1B2P", 
                        "SERPINE1", "SPP1", "STAT3", "STAT5A", "TGFB", "TGFB1", "TNFA", "TNF", 
                        "TNF-B", "TNFB", "TNFSF10", "VEGFA", "LET7", "MIR7", "MIR7-5P", "MIR15A", 
                        "MIR15B", "MIR17-92", "MIR21", "MIR23B", "MIR24", "MIR29", "MIR29A", "MIR31", 
                        "MIR33B", "MIR34B", "MIR55", "MIR75P", "MIR99A", "MIR100", "MIR125", "MIR126", 
                        "MIR127", "MIR127-3P", "MIR1273A", "MIR1273E", "MIR1273F", "MIR1273G-3P", 
                        "MIR1285-3P", "MIR1293", "MIR132", "MIR133", "MIR133A", "MIR137", "MIR138", 
                        "MIR140", "MIR145", "MIR146A", "MIR148", "MIR155", "MIR181", "MIR181B", 
                        "MIR182", "MIR183", "MIR184", "MIR185", "MIR195", "MIR200", "MIR205", "MIR206", 
                        "MIR211-3P", "MIR223", "MIR3131", "MIR3150B-3P", "MIR345", "MIR363", 
                        "MIR373-5P", "MIR375", "MIR423", "MIR4278", "MIR4430", "MIR4450", "MIR4487", 
                        "MIR4498", "MIR455-3P", "MIR455-5P", "MIR5096")
# Check immune genes
present_inflammatory_genes <- intersect(inflammatory_genes, rownames(seurat_obj))
cat("Immune genes present in the dataset:\n")
print(present_inflammatory_genes)
# Calculate immune module score
seurat_obj <- AddModuleScore(seurat_obj, features = list(inflammatory_genes), name = "inflammatory_Score")
# Visualize the module scores using a violin plot
VlnPlot(seurat_obj, features = c("inflammatory_Score1"), pt.size = 0.1, ncol = 2)
#########################
# Feature plot for the inflammatory module score
FeaturePlot(seurat_obj, 
            features = c("inflammatory_Score1"), 
            reduction = "umap", 
            pt.size = 1.5, 
            cols = c("lightgray", "blue"))
##############
# Feature plot with customized expression range
FeaturePlot(seurat_obj, 
            features = c("inflammatory_Score1"), 
            reduction = "umap", 
            pt.size = 1.5, 
            cols = c("lightgray", "blue"), 
            min.cutoff = -3, 
            max.cutoff = 1)
######################################################78
# Subset the cells belonging to cluster 1
cluster1_subset <- subset(seurat_obj, idents = 1)

# Verify the subset
DimPlot(cluster1_subset, reduction = "umap", pt.size = 2, label = TRUE)
#######################################################78
# Normalize the data for the subset
cluster1_subset <- NormalizeData(cluster1_subset, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

# Identify highly variable features
cluster1_subset <- FindVariableFeatures(cluster1_subset, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)

# Scale the data
all_genes_cluster1 <- rownames(cluster1_subset)
cluster1_subset <- ScaleData(cluster1_subset, features = all_genes_cluster1)

# Perform PCA
cluster1_subset <- RunPCA(cluster1_subset, features = VariableFeatures(object = cluster1_subset))

# Plot the elbow plot for subset (optional step)
ElbowPlot(cluster1_subset, ndims = 50)

# Recluster the subset
cluster1_subset <- FindNeighbors(cluster1_subset, dims = 1:10)
cluster1_subset <- FindClusters(cluster1_subset, resolution = 0.5)

# Run UMAP for visualization of subclusters
cluster1_subset <- RunUMAP(cluster1_subset, dims = 1:10)

# Visualize the subclusters
DimPlot(cluster1_subset, reduction = "umap", label = TRUE, pt.size = 2)
#################
# Feature plot with customized expression range
FeaturePlot(cluster1_subset, 
            features = c("inflammatory_Score1"), 
            reduction = "umap", 
            pt.size = 1.5, 
            cols = c("lightgray", "blue"), 
            min.cutoff = -3, 
            max.cutoff = 1)
###########33333333333333333333333333333333
# Update subcluster identities to make them unique (e.g., "Cluster1_Sub1", "Cluster1_Sub2", etc.)
cluster1_subset <- RenameIdents(cluster1_subset, `0` = "Cluster1_Sub1", `1` = "Cluster1_Sub2", `2` = "Cluster1_Sub3")

# Extract the subcluster identities
subcluster_labels <- Idents(cluster1_subset)

# Map the subcluster labels back into the original Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = subcluster_labels, col.name = "Subcluster")

# Update the Idents of the original Seurat object to use subclusters for cluster 1
Idents(seurat_obj, cells = colnames(cluster1_subset)) <- subcluster_labels

# Visualize the updated UMAP with new subclusters
DimPlot(seurat_obj, reduction = "umap", group.by = "Subcluster", label = TRUE, pt.size = 2)
########################
# Extract original cluster identities
original_clusters <- Idents(seurat_obj)

# Replace the identities of cluster 1 cells with their subcluster labels
original_clusters[names(subcluster_labels)] <- subcluster_labels

# Update the identities in the Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = original_clusters, col.name = "UpdatedClusters")

# Set the identities back to the updated clusters
Idents(seurat_obj) <- seurat_obj$UpdatedClusters

# Visualize the updated UMAP
DimPlot(seurat_obj, reduction = "umap", group.by = "UpdatedClusters", label = TRUE, pt.size = 2)
######################
# Enhanced Feature Plot
FeaturePlot(seurat_obj, 
            features = c("inflammatory_Score1"), 
            reduction = "umap", 
            pt.size = 2, # Slightly larger point size for better visibility
            cols = c("lightgray", "blue"), 
            min.cutoff = -3, 
            max.cutoff = 1) + 
  ggtitle("Expression of Inflammatory Score in UMAP") + # Add a plot title
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Customize title appearance
        axis.text = element_text(size = 12), # Customize axis text
        axis.title = element_text(size = 14)) # Customize axis title size
########################part1 
######################################
VlnPlot(seurat_obj, features = c("inflammatory_Score1"), group.by = "UpdatedClusters")
###############
library(dplyr)

cluster_avg_scores <- seurat_obj@meta.data %>%
  group_by(UpdatedClusters) %>%
  summarize(mean_inflammatory_score = mean(inflammatory_Score1, na.rm = TRUE))

print(cluster_avg_scores)
########################
library(dplyr)

# Fetch the module scores and cluster identities
module_scores <- FetchData(seurat_obj, vars = c("inflammatory_Score1", "UpdatedClusters"))

# Group data by clusters and calculate the average score
average_scores <- module_scores %>%
  group_by(UpdatedClusters) %>%
  summarise(AverageScore = mean(inflammatory_Score1, na.rm = TRUE))

# Print the result
print(average_scores)
##############################marking based on the cell level 
unique(seurat_obj$UpdatedClusters)
inflammatory_labels <- c("Non-inflammatory", "Inflammatory", "Non-inflammatory", 
                         "Inflammatory", "Inflammatory", "Non-inflammatory", "Non-inflammatory")
names(inflammatory_labels) <- levels(seurat_obj$UpdatedClusters)
######################
library(dplyr)

# Create a dataframe with UpdatedClusters and corresponding inflammatory labels
cluster_labels_df <- data.frame(UpdatedClusters = levels(seurat_obj$UpdatedClusters),
                                inflammatory_status = inflammatory_labels)

# Merge the inflammatory labels into Seurat metadata using UpdatedClusters
merged_metadata <- seurat_obj@meta.data %>%
  left_join(cluster_labels_df, by = c("UpdatedClusters" = "UpdatedClusters"))

# Replace Seurat metadata with the merged version
seurat_obj@meta.data <- merged_metadata
###
head(seurat_obj@meta.data)
####################
head(Embeddings(seurat_obj, reduction = "umap"))
umap_coords <- Embeddings(seurat_obj, reduction = "umap")
seurat_obj@meta.data$UMAP_1 <- umap_coords[, 1]
seurat_obj@meta.data$UMAP_2 <- umap_coords[, 2]
#########
umap_data <- FetchData(seurat_obj, vars = c("UMAP_1", "UMAP_2", "inflammatory_status"))
head(umap_data)
###############
library(dplyr)

# Calculate label positions (centroids)
label_positions <- umap_data %>%
  group_by(inflammatory_status) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# Plot UMAP with labels
DimPlot(seurat_obj, reduction = "umap", group.by = "inflammatory_status", pt.size = 2) +
  geom_text(data = label_positions, aes(x = UMAP_1, y = UMAP_2, label = inflammatory_status), size = 5, color = "black")
#####################
rownames(seurat_obj@meta.data) <- rownames(Embeddings(seurat_obj, reduction = "umap"))
all(rownames(Embeddings(seurat_obj, reduction = "umap")) == rownames(seurat_obj@meta.data))
DimPlot(seurat_obj, reduction = "umap", group.by = "inflammatory_status", label = TRUE, pt.size = 2)
###########Inflammatory: Cluster1_Sub2 (0.291), 0 (0.269), 2 (0.265)

#######Non-inflammatory: Cluster1_Sub1 (0.111), Cluster1_Sub3 (0.0685), 3 (0.0873), 4 (0.105)
#marking based on threshold 
cluster_avg_scores <- seurat_obj@meta.data %>%
  group_by(UpdatedClusters) %>%
  summarize(avg_inflammatory_score = mean(inflammatory_Score1, na.rm = TRUE))
cluster_avg_scores <- cluster_avg_scores %>%
  mutate(cluster_status = ifelse(avg_inflammatory_score > 0, "Inflammatory", "Non-inflammatory"))
seurat_obj$ClusterInflammatoryStatus <- cluster_avg_scores$cluster_status[match(seurat_obj$UpdatedClusters, cluster_avg_scores$UpdatedClusters)]
DimPlot(seurat_obj, reduction = "umap", group.by = "ClusterInflammatoryStatus", label = TRUE, pt.size = 2)
############
# Inflammatory-356  Non-inflammatory-210
cluster_labels <- data.frame(
  UpdatedClusters = c("Cluster1_Sub2", "0", "2", "Cluster1_Sub1", "Cluster1_Sub3", "3", "4"),
  cluster_status = c("Inflammatory", "Inflammatory", "Inflammatory", 
                     "Non-inflammatory", "Non-inflammatory", "Non-inflammatory", "Non-inflammatory")
)
seurat_obj$ClusterInflammatoryStatus <- cluster_labels$cluster_status[match(seurat_obj$UpdatedClusters, cluster_labels$UpdatedClusters)]
table(seurat_obj$ClusterInflammatoryStatus)
DimPlot(seurat_obj, reduction = "umap", group.by = "ClusterInflammatoryStatus", label = TRUE, pt.size = 2)
################
seurat_obj$ClusterInflammatoryStatus <- factor(
  seurat_obj$ClusterInflammatoryStatus, 
  levels = c("Inflammatory", "Non-inflammatory"),
  labels = c("Inflammatory (Cluster1_Sub2, 0, 2)", 
             "Non-inflammatory (Cluster1_Sub1, Cluster1_Sub3, 3, 4)")
)
DimPlot(seurat_obj, reduction = "umap", group.by = "ClusterInflammatoryStatus", label = FALSE, pt.size = 2) +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
##############################
umap_data <- FetchData(seurat_obj, vars = c("UMAP_1", "UMAP_2", "UpdatedClusters"))
library(dplyr)

label_positions <- umap_data %>%
  group_by(UpdatedClusters) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
DimPlot(seurat_obj, reduction = "umap", group.by = "ClusterInflammatoryStatus", label = FALSE, pt.size = 2) +
  geom_text(data = label_positions, aes(x = UMAP_1, y = UMAP_2, label = UpdatedClusters), size = 5, color = "black")
##################################################################
# Extract expression matrix
expr_matrix <- GetAssayData(seurat_obj, layer = "data")
# Extract inflammatory status metadata
inflammatory_status <- seurat_obj$ClusterInflammatoryStatus
dim(expr_matrix)  # Check dimensions of the matrix
head(rownames(expr_matrix))  # Inspect some gene names
head(colnames(expr_matrix))  # Inspect some cell barcodes
##################
library(limma)
######################################
levels(inflammatory_status) <- c("Inflammatory", "Non_inflammatory")
design <- model.matrix(~ 0 + inflammatory_status)
colnames(design) <- levels(inflammatory_status)
contrast <- makeContrasts(Inflammatory - Non_inflammatory, levels = design)
fit <- lmFit(expr_matrix, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Extract DEGs
degs <- topTable(fit, adjust = "BH", sort.by = "P", number = Inf)
head(degs)
###########
##################################
# Optional: Filter significant DEGs
# Filter significant DEGs
significant_degs <- degs[degs$adj.P.Val < 0.05, ]  # Adjusted p-value cutoff

# Assign significance categories
significant_degs$significance <- ifelse(significant_degs$adj.P.Val < 0.05 & abs(significant_degs$logFC) > 0.03,
                                        ifelse(significant_degs$logFC > 0, "Upregulated", "Downregulated"),
                                        "Not significant")

library(ggplot2)

# Identify top 5 upregulated and downregulated genes
top_upregulated <- head(significant_degs[order(-significant_degs$logFC), ], 5)
top_downregulated <- head(significant_degs[order(significant_degs$logFC), ], 5)

# Combine the top genes
top_genes <- rbind(top_upregulated, top_downregulated)

# Create a new column for labeling the top genes
significant_degs$label <- ifelse(rownames(significant_degs) %in% rownames(top_genes), rownames(top_genes), "")

# Generate volcano plot
volcano_plot <- ggplot(significant_degs, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text(aes(label = label), hjust = 0.5, vjust = -0.5, size = 3, fontface = "bold") +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log Fold Change (logFC)", y = "-log10(Adjusted P-Value)") +
  theme(plot.title = element_text(hjust = 0.5))+xlim(-2, 2) + ylim(0, 70)

# Display the plot
print(volcano_plot)
####################################
###############################

# Mark genes starting with "MIR" as green and add labels for top miRNAs
significant_degs$color <- ifelse(grepl("^MIR", rownames(significant_degs)), "green", 
                                 ifelse(significant_degs$significance == "Upregulated", "red",
                                        ifelse(significant_degs$significance == "Downregulated", "blue", "gray")))

# Identify top MIR genes for upregulated and downregulated categories
top_mir_up <- head(significant_degs[grepl("^MIR", rownames(significant_degs)) & significant_degs$significance == "Upregulated", ], 5)
top_mir_down <- head(significant_degs[grepl("^MIR", rownames(significant_degs)) & significant_degs$significance == "Downregulated", ], 5)

# Add labels for top miRNAs and other top DEGs
significant_degs$label <- ifelse(rownames(significant_degs) %in% rownames(rbind(top_mir_up, top_mir_down)), 
                                 rownames(significant_degs),
                                 ifelse(rownames(significant_degs) %in% rownames(top_genes), rownames(top_genes), ""))

# Adjust label positioning properly (downregulated labels to the left, upregulated labels to the right)
volcano_plot <- ggplot(significant_degs, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label, color = ifelse(label %in% rownames(rbind(top_mir_up, top_mir_down)), "black", NA)),
                  hjust = ifelse(significant_degs$logFC > 0,1.5,0), # Fixing the alignment for both sides
                  size = 2, fontface = "bold", max.overlaps = Inf) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray", "green" = "green", "black" = "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log Fold Change (logFC)", y = "-log10(Adjusted P-Value)") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(volcano_plot)
###############################
# Filter for miRNAs (genes starting with "MIR")
mirna_genes <- significant_degs[grepl("^MIR", rownames(significant_degs)), ]

# Separate upregulated and downregulated miRNAs
upregulated_mirna <- mirna_genes[mirna_genes$significance == "Upregulated", ]
downregulated_mirna <- mirna_genes[mirna_genes$significance == "Downregulated", ]

# Display results
cat("Upregulated miRNAs:\n")
print(upregulated_mirna)

cat("Downregulated miRNAs:\n")
print(downregulated_mirna)
###################################
library(biomaRt)
# Connect to Ensembl database for human genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Get lncRNA gene annotations
lncRNA_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                     filters = "biotype",
                     values = "lncRNA",
                     mart = ensembl)
# Check for lncRNAs in significant DEGs
significant_degs$lncRNA <- ifelse(significant_degs$label %in% lncRNA_data$external_gene_name, TRUE, FALSE)

# Count the number of lncRNAs in significant DEGs
lncRNA_count <- sum(significant_degs$lncRNA)
cat("Number of significant DEGs that are lncRNAs:\n", lncRNA_count)
#########################
library(ggplot2)
library(ggrepel)

# Mark significant lncRNAs in violet and keep other significant DEGs colored appropriately
significant_degs$color <- ifelse(significant_degs$lncRNA & significant_degs$significance == "Upregulated", "purple",
                                 ifelse(significant_degs$lncRNA & significant_degs$significance == "Downregulated", "purple",
                                        ifelse(significant_degs$significance == "Upregulated", "red",
                                               ifelse(significant_degs$significance == "Downregulated", "blue", "gray"))))

# Identify top 5 significant upregulated and downregulated lncRNAs
top_lnc_up <- head(significant_degs[significant_degs$lncRNA & significant_degs$significance == "Upregulated", ], 5)
top_lnc_down <- head(significant_degs[significant_degs$lncRNA & significant_degs$significance == "Downregulated", ], 5)

# Identify top 5 significant upregulated and downregulated normal DEGs (non-lncRNAs)
top_normal_up <- head(significant_degs[!significant_degs$lncRNA & significant_degs$significance == "Upregulated", ], 5)
top_normal_down <- head(significant_degs[!significant_degs$lncRNA & significant_degs$significance == "Downregulated", ], 5)

# Add labels for top lncRNAs and top normal DEGs
significant_degs$label <- ifelse(rownames(significant_degs) %in% rownames(rbind(top_lnc_up, top_lnc_down, top_normal_up, top_normal_down)),
                                 rownames(significant_degs), "")

# Generate the volcano plot
volcano_plot <- ggplot(significant_degs, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label),
                  color = "black",  # Black labels for top genes (lncRNAs and normal DEGs)
                  hjust = ifelse(significant_degs$logFC > 0, 0.5, 0.5),  # Adjust label positions
                  size = 3, fontface = "bold", max.overlaps = Inf) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray", "purple" = "purple"),
                     labels = c("red" = "Upregulated",
                                "blue" = "Downregulated",
                                "gray" = "Unsignificant",
                                "purple" = "lncRNAs")) +
  theme_minimal() +
  labs(title = "Volcano Plot Highlighting Significant lncRNAs and DEGs",
       x = "Log Fold Change (logFC)",
       y = "-log10(Adjusted P-Value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-2, 2) +
  ylim(0, 70)

# Display the plot
print(volcano_plot)
################
# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Annotation package for humans

# Step 1: Filter DEGs
filtered_degs <- significant_degs[significant_degs$adj.P.Val < 0.05 & abs(significant_degs$logFC) > 0.5, ]
cat("Number of filtered DEGs:\n")
print(dim(filtered_degs))

# Step 2: Prepare gene list (use gene symbols or Ensembl IDs depending on your dataset)
gene_list <- rownames(filtered_degs)

# Step 3: Perform GO enrichment analysis
go_enrichment <- enrichGO(gene         = gene_list,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",  # Change to "ENSEMBL" if using Ensembl IDs
                          ont          = "BP",     # Biological Process ontology
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# Step 4: Visualize the results
# Display results
head(go_enrichment)

# Create a bar plot of enriched GO terms
barplot(go_enrichment, showCategory = 10, title = "GO Enrichment Analysis",font.size = 8)
#####################

##########################3
# Create a dot plot for visualization
dotplot(go_enrichment, showCategory = 10, title = "GO Enrichment Analysis (all Degs)",font.size = 8)
################################
# Step 1: Filter lncRNAs from significant DEGs
lncrna_degs <- significant_degs[significant_degs$lncRNA & significant_degs$adj.P.Val < 0.05 & abs(significant_degs$logFC) > 0.5, ]
cat("Number of lncRNAs:\n")
print(dim(lncrna_degs))

# Step 2: Prepare lncRNA gene list
lncrna_genes <- rownames(lncrna_degs)

# Step 3: Run GO enrichment analysis for lncRNAs
lncrna_go_enrichment <- enrichGO(gene         = lncrna_genes,
                                 OrgDb        = org.Hs.eg.db,
                                 keyType      = "SYMBOL",  # Use "ENSEMBL" if required
                                 ont          = "BP",     # Biological Process ontology
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.2)

# Step 4: Visualize lncRNA-specific enrichment results
barplot(lncrna_go_enrichment, showCategory = 10, title = "GO Enrichment Analysis for lncRNAs",font.size = 8)
dotplot(lncrna_go_enrichment, showCategory = 10, title = "GO Enrichment Analysis for lncRNAs",font.size = 8)
#############################
# Load necessary libraries
library(pheatmap)

# Step 1: Subset top 20 DEGs
# Assuming `significant_degs` is your dataset and it contains expression values in columns
top_20_degs <- head(significant_degs[order(abs(significant_degs$logFC), decreasing = TRUE), ], 20)

# Step 2: Extract expression data
# Replace 'expression_matrix' with the actual name of the matrix containing your gene expression data
# Make sure rownames of the matrix correspond to the gene names in `top_20_degs`
expression_data <- expr_matrix[rownames(top_20_degs), ]

# Step 3: Generate Heatmap
pheatmap(expression_data,
         scale = "row",                # Scale rows to show relative expression
         cluster_rows = TRUE,          # Cluster genes
         cluster_cols = TRUE,          # Cluster samples
         show_rownames = TRUE,         # Show gene names
         show_colnames = FALSE,         # Show sample names
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
         main = "Heatmap of Top 20 DEGs")

# Optionally save the heatmap to a file
# png("heatmap_top_20_degs.png", width = 800, height = 600)
# print(pheatmap(expression_data, ...))  # Same arguments as above
# dev.off()
#################################
# Step 1: Create annotation_col using immune activation scores
annotation_col <- data.frame(InflammatoryScore = module_scores$inflammatory_Score1)
rownames(annotation_col) <- rownames(module_scores)

# Step 2: Sort samples in both expression_data and annotation_col
sorted_samples <- colnames(expression_data)[order(annotation_col$InflammatoryScore, decreasing = TRUE)]
expression_data <- expression_data[, sorted_samples]
annotation_col <- annotation_col[sorted_samples, , drop = FALSE]

# Step 3: Generate the heatmap
pheatmap(expression_data,
         fontsize_row = 8,
         scale = "row",
         cluster_rows = TRUE,          # Cluster genes
         cluster_cols = FALSE,         # Do not cluster samples (sorted by immune scores)
         annotation_col = annotation_col,  # Add immune evasion scores as annotations
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Top 20 DEGs (Sorted by Immune Evasion Scores)")
##################