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
# Subset the cells belonging to cluster 1
cluster1_subset <- subset(seurat_obj, idents = 1)

# Verify the subset
DimPlot(cluster1_subset, reduction = "umap", pt.size = 2, label = TRUE)
###
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
########################
immune_evasion_genes <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CD244", "CD160", "IRF4", "BATF", 
                             "NFATC1", "FOXP3", "IL2RA", "TNFRSF18", "IKZF4", "CD80", "CD86", "CD274", 
                             "PDCD1LG2", "LGALS9")

# Check immune genes
present_immune_evasion_genes <- intersect(immune_evasion_genes, rownames(seurat_obj))
cat("Immune evasion genes present in the dataset:\n")
print(present_immune_evasion_genes)
dim(seurat_obj)
####################################################################
# Calculate immune evasion module score
seurat_obj <- AddModuleScore(seurat_obj, features = list(immune_evasion_genes), name = "Immune_evasion_Score")


# Visualize the module scores using a violin plot
VlnPlot(seurat_obj, features = c("Immune_evasion_Score1"), pt.size = 0, ncol = 2)
#########################
library(dplyr)

# Fetch the module scores and cluster identities
module_scores <- FetchData(seurat_obj, vars = c("Immune_evasion_Score1", "UpdatedClusters"))

# Group data by clusters and calculate the average score
average_scores <- module_scores %>%
  group_by(UpdatedClusters) %>%
  summarise(AverageScore = mean(Immune_evasion_Score1, na.rm = TRUE))

# Print the result
print(average_scores)
############################################
# Extract normalized gene expression data using GetAssayData
expression_data <- as.matrix(GetAssayData(seurat_obj, slot = "data"))

# Verify the dimensions of the extracted data
cat("Expression data dimensions:\n")
print(dim(expression_data))
##############

# Step 1: Extract normalized gene expression data and updated cluster information
expression_data <- as.matrix(GetAssayData(seurat_obj, slot = "data"))  # Gene expression matrix
updated_clusters <- seurat_obj@meta.data$UpdatedClusters              # Updated cluster labels

# Step 2: Define groups for comparison using updated clusters
group <- ifelse(updated_clusters %in% c("Cluster1_Sub2", "Cluster1_Sub1", "2"), "GroupA", "GroupB")

# Step 3: Create a design matrix
library(limma)
group <- factor(group, levels = c("GroupA", "GroupB"))                 # Ensure group is a factor
design <- model.matrix(~ 0 + group)
colnames(design) <- c("GroupA", "GroupB")                              # Define column names

# Verify the design matrix
print(design)


# Step 4: Fit a linear model
fit <- lmFit(expression_data, design)

# Step 5: Set up contrast for comparison
contrast <- makeContrasts(GroupA - GroupB, levels = design)

# Step 6: Apply the contrast and calculate DEGs
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Step 7: Extract results
deg_results <- topTable(fit, number = Inf, sort.by = "p")
cat("Top DEGs:\n")
print(head(deg_results))  # View the top DEGs
dim(deg_results)
# Optional: Filter significant DEGs
# Filter significant DEGs
significant_degs <- deg_results[deg_results$adj.P.Val < 0.05, ]  # Adjusted p-value cutoff

# Assign significance categories
significant_degs$significance <- ifelse(significant_degs$adj.P.Val < 0.05 & abs(significant_degs$logFC) > 0.03,
                                        ifelse(significant_degs$logFC > 0, "Upregulated", "Downregulated"),
                                        "Not significant")

# Display significant DEGs
cat("Significant DEGs with significance category:\n")
print(head(significant_degs))


# Step 8: Export significant DEGs to a CSV file (optional)
write.csv(significant_degs, file = "Significant_DEGs_GroupA_vs_GroupB.csv", row.names = TRUE)
cat("Significant DEGs saved to 'Significant_DEGs_GroupA_vs_GroupB.csv'\n")
#########################################
# Check if row names contain Ensembl IDs
ensembl_ids <- rownames(deg_results)

# Count the total number of Ensembl IDs
num_ensembl_ids <- sum(grepl("^ENS", ensembl_ids))  # Ensembl IDs typically start with "ENS"

# Print the number of Ensembl IDs
cat("Number of Ensembl IDs in DEGs data:", num_ensembl_ids, "\n")
dim(deg_results)
##################################
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
library(ggplot2)
library(ggrepel)

# Mark genes starting with "MIR" as green and add labels for top miRNAs
significant_degs$color <- ifelse(grepl("^MIR", rownames(significant_degs)), "green", 
                                 ifelse(significant_degs$significance == "Upregulated", "red",
                                        ifelse(significant_degs$significance == "Downregulated", "blue", "gray")))

# Identify top MIR genes for upregulated and downregulated categories
top_mir_up <- head(significant_degs[grepl("^MIR", rownames(significant_degs)) & significant_degs$significance == "Upregulated", ], 5)
top_mir_down <- head(significant_degs[grepl("^MIR", rownames(significant_degs)) & significant_degs$significance == "Downregulated", ],5)

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
########################################
# Filter for miRNAs (genes starting with "MIR")
mirna_genes <- significant_degs[grepl("^MIR", rownames(significant_degs)), ]
# Check the number of miRNAs found
num_mirna <- nrow(mirna_genes)

# Print the result
cat("Number of miRNAs found:\n", num_mirna)
# Separate upregulated and downregulated miRNAs
upregulated_mirna <- mirna_genes[mirna_genes$significance == "Upregulated", ]
downregulated_mirna <- mirna_genes[mirna_genes$significance == "Downregulated", ]

# Display results
cat("Upregulated miRNAs:\n")
print(upregulated_mirna)

cat("Downregulated miRNAs:\n")
print(downregulated_mirna)
#################################################don't use 
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
significant_degs$color <- ifelse(significant_degs$lncRNA & significant_degs$significance == "Upregulated", "violet",
                                 ifelse(significant_degs$lncRNA & significant_degs$significance == "Downregulated", "violet",
                                        ifelse(significant_degs$significance == "Upregulated", "red",
                                               ifelse(significant_degs$significance == "Downregulated", "blue", "gray"))))

# Identify top significant upregulated and downregulated lncRNAs
top_lnc_up <- head(significant_degs[significant_degs$lncRNA & significant_degs$significance == "Upregulated", ], 5)
top_lnc_down <- head(significant_degs[significant_degs$lncRNA & significant_degs$significance == "Downregulated", ], 5)

# Add labels for top lncRNAs
significant_degs$label <- ifelse(rownames(significant_degs) %in% rownames(rbind(top_lnc_up, top_lnc_down)),
                                 rownames(significant_degs), "")

# Generate the volcano plot
volcano_plot <- ggplot(significant_degs, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label),
                  color = "black",  # Black labels for top lncRNAs
                  hjust = ifelse(significant_degs$logFC > 0, 1.5, 0),  # Adjust label positions
                  size = 3, fontface = "bold", max.overlaps = Inf) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray", "violet" = "violet")) +
  theme_minimal() +
  labs(title = "Volcano Plot Highlighting Significant lncRNAs",
       x = "Log Fold Change (logFC)",
       y = "-log10(Adjusted P-Value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-2, 2) +
  ylim(0, 70)

# Display the plot
print(volcano_plot)

#####################################
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
                  hjust = ifelse(significant_degs$logFC > 0, 1.5, 0),  # Adjust label positions
                  size = 3, fontface = "bold", max.overlaps = Inf) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray", "purple" = "purple")) +
  theme_minimal() +
  labs(title = "Volcano Plot Highlighting Significant lncRNAs and DEGs",
       x = "Log Fold Change (logFC)",
       y = "-log10(Adjusted P-Value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-2, 2) +
  ylim(0, 70)

# Display the plot
print(volcano_plot)

##############################################$$$$$
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
                  hjust = ifelse(significant_degs$logFC > 0, 1.2, 0.5),  # Adjust label positions
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

################################################
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
######################
library(clusterProfiler)
library(ggplot2)

# Extract the enrichment results data
go_data <- as.data.frame(go_enrichment)

# Generate the bar plot with annotations for p-values
barplot <- ggplot(go_data[1:10, ], aes(x = reorder(Description, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = ifelse(p.adjust < 0.001, "p < 0.001", paste0("p = ", round(p.adjust, 3)))),
            vjust = -0.5, size = 0.1) +  # Annotate with simplified p-values
  theme_minimal() +
  labs(title = "GO Enrichment Analysis",
       x = "GO Terms",
       y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(barplot)

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
expression_data <- expression_data[rownames(top_20_degs), ]

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
# Step 1: Create annotation_col using immune evasion scores
annotation_col <- data.frame(ImmuneEvasionScore = module_scores$Immune_evasion_Score1)
rownames(annotation_col) <- rownames(module_scores)

# Step 2: Sort samples in both expression_data and annotation_col
sorted_samples <- colnames(expression_data)[order(annotation_col$ImmuneEvasionScore, decreasing = TRUE)]
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
