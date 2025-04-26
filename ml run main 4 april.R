working_data <- readRDS("working_data.rds")
str(working_data)
###
# Load the required library
library(glmnet)

# Feature matrix: Convert gene expression data to numeric and ensure samples are columns
x <- t(as.matrix(working_data[-which(rownames(working_data) == "Immune_evasion_Label"), ]))
x <- as.numeric(x) # Convert the matrix to numeric
x <- matrix(x, nrow = ncol(working_data), ncol = nrow(working_data) - 1) # Reshape dimensions (samples as rows)

# Restore row and column names (ensure gene names and sample identifiers are intact)
colnames(x) <- rownames(working_data)[-which(rownames(working_data) == "Immune_evasion_Label")] # Gene names
rownames(x) <- colnames(working_data) # Sample identifiers

# Response vector (immune evasion labels corresponding to samples)
y <- as.numeric(working_data["Immune_evasion_Label", ]) # Dependent variable (evasive or non-evasive)

# Ensure dimensions match between feature matrix and response vector
cat("Dimensions of x (samples as rows, genes as columns):", dim(x), "\n")
cat("Length of y:", length(y), "\n")

# Check for missing values in the matrix and response vector
if (sum(is.na(x)) > 0 || sum(is.na(y)) > 0) {
  stop("Missing values detected in x or y. Please clean the data before proceeding.")
}

# Standardize the feature matrix (scaling gene expression values)
x <- scale(x)

# Perform LASSO regression with cross-validation
lasso_model <- cv.glmnet(x, y, alpha = 1, family = "binomial")

# Extract non-zero coefficients (selected features)
selected_features <- which(coef(lasso_model, s = "lambda.min") != 0)

# Map selected features back to gene names (excluding intercept)
selected_genes <- colnames(x)[selected_features[-1]] # Exclude the intercept
print(selected_genes)

# Save selected genes for future use
saveRDS(selected_genes, file = "selected_genes_lasso.rds")

# Print confirmation
cat("Selected genes using LASSO regression:", length(selected_genes), "genes\n")
#####################################3
selected_genes <- readRDS(file = "selected_genes_lasso.rds")

# Load necessary libraries
library(caret)
library(glmnet)
library(pROC)
library(ggplot2)

# --------------------------
# Assume "model_data" has been built as follows:
#   model_data: a data frame with selected gene features (columns) and the response variable "Immune_evasion_Label"
# In our case, we have:
#    x: samples (rows) x genes (columns) from previous steps (and standardized)
#    y: response vector (0/1) converted to factor ("NonEvasive", "Evasive")
# and we created model_data from the selected_genes:
model_data <- as.data.frame(x[, selected_genes, drop = FALSE])
model_data$Immune_evasion_Label <- factor(y, levels = c(0, 1), labels = c("NonEvasive", "Evasive"))
cat("Dimensions of model_data (samples x features):", dim(model_data), "\n")

# --------------------------
# 1. Partition the Dataset (80% train, 20% test)
set.seed(123)
trainIndex <- createDataPartition(model_data$Immune_evasion_Label, p = 0.8, list = FALSE)
train_data <- model_data[trainIndex, ]
test_data  <- model_data[-trainIndex, ]
cat("Training set dimensions:", dim(train_data), "\n")
cat("Testing set dimensions:", dim(test_data), "\n")

# --------------------------
# 2. Set up 10-fold Cross-Validation (optimizing ROC)
train_control <- trainControl(method = "cv",
                              number = 10,
                              classProbs = TRUE,
                              summaryFunction = twoClassSummary,
                              savePredictions = TRUE)

# --------------------------
# 3. Build ML Models

# Logistic Regression
set.seed(123)
logistic_model <- train(Immune_evasion_Label ~ .,
                        data = train_data,
                        method = "glm",
                        family = "binomial",
                        trControl = train_control,
                        metric = "ROC")
cat("Logistic Regression Model Summary:\n")
print(logistic_model)

# Support Vector Machine with Radial Kernel
set.seed(123)
svm_model <- train(Immune_evasion_Label ~ .,
                   data = train_data,
                   method = "svmRadial",
                   preProcess = c("center", "scale"),
                   trControl = train_control,
                   metric = "ROC",
                   tuneLength = 5)
cat("SVM Model Summary:\n")
print(svm_model)

# Random Forest
set.seed(123)
rf_model <- train(Immune_evasion_Label ~ .,
                  data = train_data,
                  method = "rf",
                  trControl = train_control,
                  metric = "ROC",
                  tuneLength = 5)
cat("Random Forest Model Summary:\n")
print(rf_model)

# --------------------------
# 4. Evaluate Models on the Test Set

# For Logistic Regression:
logistic_preds_prob <- predict(logistic_model, newdata = test_data, type = "prob")
logistic_preds_class <- predict(logistic_model, newdata = test_data)
cm_logistic <- confusionMatrix(logistic_preds_class, test_data$Immune_evasion_Label)
roc_logistic <- roc(response = test_data$Immune_evasion_Label,
                    predictor = logistic_preds_prob$Evasive,
                    levels = c("NonEvasive","Evasive"),
                    direction = "<")
auc_logistic <- auc(roc_logistic)

# For SVM:
svm_preds_prob <- predict(svm_model, newdata = test_data, type = "prob")
svm_preds_class <- predict(svm_model, newdata = test_data)
cm_svm <- confusionMatrix(svm_preds_class, test_data$Immune_evasion_Label)
roc_svm <- roc(response = test_data$Immune_evasion_Label,
               predictor = svm_preds_prob$Evasive,
               levels = c("NonEvasive","Evasive"),
               direction = "<")
auc_svm <- auc(roc_svm)

# For Random Forest:
rf_preds_prob <- predict(rf_model, newdata = test_data, type = "prob")
rf_preds_class <- predict(rf_model, newdata = test_data)
cm_rf <- confusionMatrix(rf_preds_class, test_data$Immune_evasion_Label)
roc_rf <- roc(response = test_data$Immune_evasion_Label,
              predictor = rf_preds_prob$Evasive,
              levels = c("NonEvasive","Evasive"),
              direction = "<")
auc_rf <- auc(roc_rf)

# --------------------------
# 5. Compare Model Performance
performance <- data.frame(
  Model = c("Logistic Regression", "SVM", "Random Forest"),
  Accuracy = c(cm_logistic$overall["Accuracy"],
               cm_svm$overall["Accuracy"],
               cm_rf$overall["Accuracy"]),
  Sensitivity = c(cm_logistic$byClass["Sensitivity"],
                  cm_svm$byClass["Sensitivity"],
                  cm_rf$byClass["Sensitivity"]),
  Specificity = c(cm_logistic$byClass["Specificity"],
                  cm_svm$byClass["Specificity"],
                  cm_rf$byClass["Specificity"]),
  AUC = c(auc_logistic, auc_svm, auc_rf)
)
print(performance)

# Print confusion matrices for each model:
cat("Confusion Matrix for Logistic Regression:\n")
print(cm_logistic)
cat("Confusion Matrix for SVM:\n")
print(cm_svm)
cat("Confusion Matrix for Random Forest:\n")
print(cm_rf)

# --------------------------
# 6. Plot ROC Curves for All Models
plot(roc_logistic, col = "red", lwd = 2, main = "ROC Curves Comparison")
plot(roc_svm, col = "blue", lwd = 2, add = TRUE)
plot(roc_rf, col = "green", lwd = 2, add = TRUE)
legend("bottomright", legend = c(paste("Logistic (AUC =", round(auc_logistic, 3), ")"),
                                 paste("SVM (AUC =", round(auc_svm, 3), ")"),
                                 paste("RF (AUC =", round(auc_rf, 3), ")")),
       col = c("red", "blue", "green"), lwd = 2)
###############################
# Check the column names and ensure they match what the model expects
print(colnames(X_train))

# Check the overall structure of the data 
str(X_train)

# See a summary of each variable to look for missing values or unusual values
summary(X_train)

  


#################################
# --- Re-fit the Model Using the x, y Interface ---
# 1. Prepare X_train (all predictor columns)
X_train <- train_data[, setdiff(names(train_data), "Immune_evasion_Label"), drop = FALSE]
colnames(X_train) <- make.names(colnames(X_train))  # Ensure sanitized names

# 2. Fit the model with x and y
library(randomForest)
rf_model <- randomForest(x = X_train,
                         y = as.factor(train_data$Immune_evasion_Label))

# --- Now proceed with SHAP analysis ---
library(fastshap)
library(ggplot2)

# In this case, since we used the x, y interface,
# the model's stored predictor names should match X_train.
expected_vars <- colnames(X_train)

# Diagnostic check:
missing_vars <- setdiff(expected_vars, colnames(X_train))
if (length(missing_vars) > 0) {
  cat("Missing variables:", missing_vars, "\n")
} else {
  cat("All expected variables are present in X_train.\n")
}

# Use the same X_train
# Define a prediction wrapper with robust imputation.
rf_pred_impute <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  complete_newdata <- as.data.frame(matrix(NA_real_, nrow = nrow(newdata), ncol = length(expected_vars)))
  colnames(complete_newdata) <- expected_vars
  
  for (var in expected_vars) {
    if (var %in% colnames(newdata)) {
      complete_newdata[[var]] <- newdata[[var]]
    } else {
      cat("Variable", var, "is missing in newdata. Imputing median.\n")
      complete_newdata[[var]] <- rep(median(X_train[[var]], na.rm = TRUE), nrow(newdata))
    }
  }
  
  for (j in seq_along(complete_newdata)) {
    if (any(is.na(complete_newdata[[j]]))) {
      complete_newdata[[j]][is.na(complete_newdata[[j]])] <- median(X_train[[j]], na.rm = TRUE)
    }
  }
  
  predict(object, newdata = complete_newdata, type = "prob")[, "Evasive"]
}

# Compute SHAP values
set.seed(123)
shap_rf <- explain(object = rf_model,  # use the re-fitted rf_model
                   X = X_train, 
                   pred_wrapper = rf_pred_impute, 
                   nsim = 50,
                   adjust = TRUE)

# Summary Plot
mean_abs_shap <- colMeans(abs(shap_rf))
df_shap <- data.frame(Feature = names(mean_abs_shap), MeanAbsSHAP = mean_abs_shap)
df_shap <- df_shap[order(df_shap$MeanAbsSHAP, decreasing = TRUE), ]

ggplot(df_shap, aes(x = reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "SHAP Feature Importance (Random Forest)",
       x = "Feature (Gene)",
       y = "Mean Absolute SHAP Value") +
  theme_minimal() +  # Note the parentheses
  theme(
    axis.title = element_text(size = 10),   # Adjust axis title size
    axis.text  = element_text(size = 8),      # Adjust axis tick text size
    plot.title = element_text(size = 12)      # Adjust plot title size
  )
################################
print(df_shap)

# Subset to the top 20 features
df_shap_top20 <- df_shap[1:20, ]

# Plot the top 20 features
ggplot(df_shap_top20, aes(x = reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "SHAP Feature Importance (Random Forest)",
       x = "Feature (Gene)",
       y = "Mean Absolute SHAP Value") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),   # Adjust axis title size
    axis.text  = element_text(size = 8),      # Adjust axis tick text size
    plot.title = element_text(size = 12)      # Adjust plot title size
  )
###############################################
library(biomaRt)

# Connect to Ensembl database for human genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get lncRNA gene annotations
lncRNA_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                     filters = "biotype",
                     values = "lncRNA",
                     mart = ensembl)

# Convert row names of df_shap and external_gene_name in lncRNA_data to uppercase for consistency
rownames(df_shap) <- toupper(rownames(df_shap))
lncRNA_data$external_gene_name <- toupper(lncRNA_data$external_gene_name)

# Check for lncRNAs in df_shap using row names
df_shap$lncRNA <- ifelse(rownames(df_shap) %in% lncRNA_data$external_gene_name, TRUE, FALSE)

# Create an object with only the lncRNA rows
lncRNA_shap <- df_shap[df_shap$lncRNA == TRUE, ]

# View the dimensions and preview the lncRNA_shap object
dim(lncRNA_shap)
head(lncRNA_shap)

# Save the lncRNA subset object for future use
saveRDS(lncRNA_shap, "lncRNA_shap.rds")

# Print confirmation
cat("lncRNA subset has been created and saved as 'lncRNA_shap.rds'\n")

#########################################
# OPTIONAL: Dependence plot for "ACTB"
feature_to_plot <- "ACTB"
plot_df <- data.frame(FeatureValue = X_train[[feature_to_plot]],
                      SHAP = shap_rf[, feature_to_plot])
ggplot(plot_df, aes(x = FeatureValue, y = SHAP)) +
  geom_point(alpha = 0.5, color = "darkorange") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(title = paste("Dependence Plot for", feature_to_plot),
       x = feature_to_plot,
       y = "SHAP Value") +
  theme_minimal()
##############
print((df_shap_top20))
######################
output_dir <- "C:\\Users\\KIIT\\OneDrive\\Documents\\peter_pan"
# List of top 5 genes
genes <- c("CD74", "CYBB", "SECTM1", "IGHM", "C1R")

# Loop through each gene and create the dependence plot
for (feature_to_plot in genes) {
  # Create the dataframe for the current gene
  plot_df <- data.frame(
    FeatureValue = X_train[[feature_to_plot]],
    SHAP = shap_rf[, feature_to_plot]
  )
  
  # Generate the ggplot for the current gene
  p <- ggplot(plot_df, aes(x = FeatureValue, y = SHAP)) +
    geom_point(alpha = 0.5, color = "darkorange") +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(
      title = paste("Dependence Plot for", feature_to_plot),
      x = feature_to_plot,
      y = "SHAP Value"
    ) +
    theme_minimal()
  
  # Print or save the plot (optional)
  print(p)
  
  # Optionally, save to a file
  # Save the plot to the specified location
  ggsave(filename = paste0(output_dir, "\\", feature_to_plot, "_dependence_plot.png"), plot = p)
  }

############################
###
library(GSVA)          # Loads the GSVA package
library(BiocParallel)  # Loads the BiocParallel package
library(ReactomePA)    # Loads the ReactomePA package
################
library(org.Hs.eg.db)

# Map gene symbols to Entrez IDs
entrez_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = df_shap$Feature,       # Use the Feature column (gene symbols)
  column = "ENTREZID",          # Convert to Entrez IDs
  keytype = "SYMBOL",           # Original key type is SYMBOL
  multiVals = "first"           # Handle duplicates by selecting the first mapping
)

# Add Entrez IDs to df_shap
df_shap$EntrezID <- entrez_ids

# Filter genes with valid Entrez IDs
df_shap_filtered <- df_shap[!is.na(df_shap$EntrezID), ]

# Create ranked gene list for GSEA
gene_list <- setNames(df_shap_filtered$MeanAbsSHAP, df_shap_filtered$EntrezID)
################
range(gene_list)
library(clusterProfiler)
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,  # Human database
  keyType = "ENTREZID",
  ont = "BP",            # Biological process
  minGSSize = 10,
  pvalueCutoff = 0.05,
  scoreType = "pos",     # Specify that all scores are positive
  verbose = TRUE
)
library(enrichplot)
# Visualize the top 10 enriched GO terms
dotplot(gsea_go, showCategory = 10) + ggtitle("GSEA GO Dotplot")
# Check the top entries to decide which gene set to plot:
head(gsea_go@result)

ridgeplot(gsea_go) +
  ggtitle("GSEA GO Ridge Plot") +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)
  )

cnetplot(gsea_go, showCategory = 5) + ggtitle("GSEA GO Category Network")
nrow(gsea_go@result)
# Convert the enrichment result to a readable format by mapping Entrez IDs to gene symbols
gsea_go_readable <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Create the cnet plot with gene names instead of numeric IDs
cnetplot(gsea_go_readable, showCategory = 5) + 
  ggtitle("GSEA GO Category Network")
################################
library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)

# Get mapping of Entrez IDs to gene symbols and gene function names
gene_info <- select(org.Hs.eg.db,
                    keys = df_shap_filtered$EntrezID,
                    columns = c("SYMBOL", "GENENAME"),
                    keytype = "ENTREZID")

# Convert enrichment result to a readable format (this maps to gene symbols)
gsea_go_readable <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Function to replace gene symbols with gene function names using 'gene_info'
replace_gene_symbols <- function(geneIDs, mapping) {
  genes <- unlist(str_split(geneIDs, "/"))
  genenames <- sapply(genes, function(g) {
    idx <- which(mapping$SYMBOL == g)
    if (length(idx) > 0) {
      return(mapping$GENENAME[idx[1]])
    } else {
      return(g)
    }
  })
  paste(genenames, collapse = "/")
}

# Update the 'geneID' field in the enrichment object's result with function names
gsea_go_readable@result$geneID <- sapply(gsea_go_readable@result$geneID, 
                                         function(x) replace_gene_symbols(x, gene_info))

# Now generate the cnet plot with gene function names in place of symbols
cnetplot(gsea_go_readable, showCategory = 5) + 
  ggtitle("GSEA GO Category Network with Gene Function Names")

##########################################
# Check if you have multiple enriched terms
if (nrow(gsea_go@result) > 1) {
  # Compute pairwise similarities
  sim <- pairwise_termsim(gsea_go)
  
  # Plot the enrichment map using the computed similarities
  p <- emapplot(sim) + ggtitle("Enrichment Map")
  print(p)
} else {
  message("Not enough enriched terms for an enrichment map.")
}
###########
library(ggrepel)
update_geom_defaults("text_repel", list(size = 2.5))

# Check if you have multiple enriched terms
if (nrow(gsea_go@result) > 1) {
  # Compute pairwise similarities
  sim <- pairwise_termsim(gsea_go)
  
  # Plot the enrichment map using the computed similarities,
  # add a specific color gradient based on 'p.adjust',
  # and reduce text sizes via theme modifications
  p <- emapplot(sim, color = "p.adjust") +
    ggtitle("Enrichment Map") +
    scale_color_gradientn(
      colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027")
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 8)
    )
  
  print(p)
} else {
  message("Not enough enriched terms for an enrichment map.")
}
