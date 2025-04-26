
# Define the file path
file_path <- "C:/Users/KIIT/Downloads/oral cancer/Tumornames.csv"

# Load the file
tumor_data <- read.csv(file_path)

# Check the first few rows of the data
dim(tumor_data)

######################
library(biomaRt)

# Connect to Ensembl database for human genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get lncRNA gene annotations
lncRNA_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                     filters = "biotype",
                     values = "lncRNA",
                     mart = ensembl)

# Convert row names of deg_results and external_gene_name in lncRNA_data to uppercase for consistency
rownames(deg_results) <- toupper(rownames(deg_results))
lncRNA_data$external_gene_name <- toupper(lncRNA_data$external_gene_name)

# Check for lncRNAs in deg_results using row names
deg_results$lncRNA <- ifelse(rownames(deg_results) %in% lncRNA_data$external_gene_name, TRUE, FALSE)

# Create an object with only the lncRNA DEGs
lncRNA_degs <- deg_results[deg_results$lncRNA == TRUE, ]

# View the dimensions and preview the lncRNA_degs object
dim(lncRNA_degs)
head(lncRNA_degs)

# Save the lncRNA DEGs object for future use
saveRDS(lncRNA_degs, "lncRNA_degs.rds")

# Print confirmation
cat("lncRNA DEGs subset has been created and saved as 'lncRNA_degs.rds'\n")

###########################################################
# Ensure consistency in row names of lncRNA_degs and the X column in tumor_data
rownames(lncRNA_degs) <- toupper(rownames(lncRNA_degs))
tumor_data$X <- toupper(tumor_data$X)

# Merge lncRNA_degs with tumor_data using the X column for gene name comparison
lncRNA_degs_ml <- merge(lncRNA_degs, tumor_data, by.x = "row.names", by.y = "X", all.x = TRUE)

# Restore row names after merging
rownames(lncRNA_degs_ml) <- lncRNA_degs_ml$Row.names
lncRNA_degs_ml$Row.names <- NULL

# Check the first few rows of the merged dataset
head(lncRNA_degs_ml)
##############
# Identify rows where gene names start with "MIR"
rows_to_remove <- grepl("^MIR", rownames(lncRNA_degs_ml))

# Remove these rows from the dataset
lncRNA_degs_ml <- lncRNA_degs_ml[!rows_to_remove, ]

# Verify the dimensions and preview the updated dataset
cat("Dimensions of lncRNA_degs_ml after removing genes starting with 'MIR':", dim(lncRNA_degs_ml), "\n")
head(lncRNA_degs_ml)
####
# Count rows where the gene name (row names) starts with "MIR"
mirna_lncRNA_count <- sum(grepl("^MIR", rownames(lncRNA_degs_ml)))

# Print the result
cat("Number of gene names in lncRNA_degs_ml starting with 'MIR':", mirna_lncRNA_count, "\n")


# Save the updated dataset
saveRDS(lncRNA_degs_ml, "lncRNA_degs_ml_updated.rds")

# Confirmation
cat("Rows with gene names starting with 'MIR' have been removed and the updated dataset is saved.\n")


# Save the ML-ready dataset
saveRDS(lncRNA_degs_ml, "lncRNA_degs_for_ml_ready.rds")

# Print confirmation and dimensions of the final dataset
cat("lncRNA DEGs dataset is now ML-ready and has been saved!\n")
cat("Dimensions of the ML-ready dataset:", dim(lncRNA_degs_ml), "\n")

###############################################################
# Filter row names in deg_results that start with "MIR"
mirna_degs <- deg_results[grepl("^MIR", rownames(deg_results)), ]

# Count the number of rows starting with "MIR"
mirna_count <- nrow(mirna_degs)
cat("Number of DEGs with names starting with 'MIR':", mirna_count, "\n")

# View the first few entries in mirna_degs for validation
head(mirna_degs)
###
# Convert row names in mirna_degs and the X column of tumor_data to uppercase for consistent matching
rownames(mirna_degs) <- toupper(rownames(mirna_degs))
tumor_data$X <- toupper(tumor_data$X)

# Subset tumor_data to match row names from mirna_degs using X column
mirna_data_merged <- tumor_data[tumor_data$X %in% rownames(mirna_degs), ]

# Merge gene expression and samples with mirna_degs
mirna_degs_ml <- merge(mirna_degs, mirna_data_merged, by.x = "row.names", by.y = "X", all.x = TRUE)

# Restore row names after merging
rownames(mirna_degs_ml) <- mirna_degs_ml$Row.names
mirna_degs_ml$Row.names <- NULL
# Count the number of DEGs where lncRNA is TRUE in mirna_degs
lncRNA_true_count <- sum(mirna_degs$lncRNA == TRUE, na.rm = TRUE)

# Print the count
cat("Number of DEGs with lncRNA marked as TRUE:", lncRNA_true_count, "\n")

# Save the ML-ready dataset
saveRDS(mirna_degs_ml, "mirna_data_for_ml_ready.rds")

# View confirmation and dataset information
cat("Mirna DEGs dataset is now ML-ready and has been saved!\n")
cat("Dimensions of the ML-ready dataset:", dim(mirna_degs_ml), "\n")
colnames(mirna_degs_ml,10)
sum(is.na(cleaned_data))
lncrna_degs
###########################
# Check for duplicate columns
duplicate_columns <- duplicated(colnames(mirna_degs_ml))
cat("Duplicate columns found:\n")
print(colnames(mirna_degs_ml)[duplicate_columns])

# Check for duplicate rows
duplicate_rows <- duplicated(mirna_degs_ml)
cat("Number of duplicate rows found:", sum(duplicate_rows), "\n")

# View the duplicate rows if needed
if (sum(duplicate_rows) > 0) {
  cat("Duplicate rows:\n")
  print(mirna_degs_ml[duplicate_rows, ])
}

##################################
# Convert all row names to uppercase for consistency
rownames(deg_results) <- toupper(rownames(deg_results))
rownames(mirna_degs) <- toupper(rownames(mirna_degs))
rownames(lncRNA_degs) <- toupper(rownames(lncRNA_degs))

# Identify rows to remove
rows_to_remove <- union(rownames(mirna_degs), rownames(lncRNA_degs))

# Filter out matching rows from deg_results
filtered_degs <- deg_results[!(rownames(deg_results) %in% rows_to_remove), ]
dim(filtered_degs)
# Save the filtered dataset
saveRDS(filtered_degs, "filtered_degs_for_analysis.rds")

# Print the dimensions of the filtered dataset
cat("Number of DEGs after filtering:", nrow(filtered_degs), "\n")
######################################
# Ensure consistency in row names of filtered_degs and the X column of tumor_data
rownames(filtered_degs) <- toupper(rownames(filtered_degs))
tumor_data$X <- toupper(tumor_data$X)

# Merge filtered_degs with tumor_data using the X column for gene name comparison
filtered_degs_ml <- merge(filtered_degs, tumor_data, by.x = "row.names", by.y = "X", all.x = TRUE)

# Restore row names after merging
rownames(filtered_degs_ml) <- filtered_degs_ml$Row.names
filtered_degs_ml$Row.names <- NULL
head(filtered_degs_ml)

# Save the ML-ready dataset
saveRDS(filtered_degs_ml, "filtered_degs_for_ml_ready.rds")

# Print confirmation and dimensions of the final dataset
cat("Filtered DEGs dataset is now ML-ready and has been saved!\n")
cat("Dimensions of the ML-ready dataset:", dim(filtered_degs_ml), "\n")
#
