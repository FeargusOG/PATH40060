library(dplyr)
library(survival)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")

library(qvalue)

prepare_patient_survivability <- function(patient_data) {
  ## Set the PATIENT_ID column as the rownames
  Y_patient <- data.frame(
    patient_id = patient_data$PATIENT_ID,
    age = patient_data$AGE,
    time = patient_data$OS_MONTHS,
    status = patient_data$OS_STATUS
  )
  rownames(Y_patient) <- Y_patient$patient_id
  
  ## Drop the patient_id column
  Y_patient$patient_id <- NULL
  ## Ensure the OS_MONTHS column is numeric and > 0
  Y_patient$time <- as.numeric(Y_patient$time)
  Y_patient <- Y_patient[Y_patient$time > 0,]
  ## Conver the OS_STATUS to numeric 1 or 0
  Y_patient$status <- ifelse(Y_patient$status == "1:DECEASED", 1, 
                             ifelse(Y_patient$status == "0:LIVING", 0, NA))
  
  return(Y_patient)
}

# Load the clinical data
patient_data <- read.delim("data/gbm_tcga_pub2013/data_clinical_patient.txt", skip = 4)
survive_data <- prepare_patient_survivability(patient_data)

# Load the mRNA data
rna_data <- read.delim("data/gbm_tcga_pub2013/data_mrna_seq_v2_rsem.txt")

# Remove duplicates: keep the row with highest total expression for each Entrez ID
rna_data_filtered <- rna_data %>%
  mutate(row_sum = rowSums(across(where(is.numeric)))) %>%
  group_by(Entrez_Gene_Id) %>%
  slice_max(order_by = row_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-row_sum)

# Now transpose to get samples as rows
rna_t <- t(rna_data_filtered[, -c(1, 2)])  # Remove Hugo_Symbol and Entrez_Gene_Id before transpose
rna_t <- as.data.frame(rna_t)
colnames(rna_t) <- rna_data_filtered$Entrez_Gene_Id

# Clean sample IDs to match survival data
cleaned_rownames <- gsub("\\.\\d+$", "", rownames(rna_t))
cleaned_rownames <- gsub("\\.", "-", cleaned_rownames)
rownames(rna_t) <- cleaned_rownames

# Ensure samples match
common_samples <- intersect(rownames(survive_data), rownames(rna_t))
surv <- survive_data[common_samples, ]
rna_t <- rna_t[common_samples, ]

# ----- Filtering Parameters -----
min_variance <- 0.6             # Minimum acceptable variance for a gene
min_expr_prop <- 0.2            # Minimum proportion of samples in which gene must be expressed (> 0)

# ----- Filtering Logic -----

# 1. Calculate variance for each gene
gene_variance <- apply(rna_t, 2, var, na.rm = TRUE)

# 2. Calculate expression frequency for each gene
#    (i.e., proportion of samples where expression > 0)
gene_expr_prop <- colMeans(rna_t > 0)

# 3. Identify genes that pass both filters
keep_genes <- (gene_variance > min_variance) & (gene_expr_prop > min_expr_prop)

# 4. Subset the RNA data to only include these genes
rna_t_filtered <- rna_t[, keep_genes]

# z-score scale
rna_scaled <- as.data.frame(scale(rna_t_filtered))

# Run Cox model for each gene
results <- lapply(colnames(rna_scaled), function(gene) {
  df <- cbind(surv, gene_expr = rna_scaled[[gene]])
  model <- coxph(Surv(time, status) ~ gene_expr + age, data = df)
  summary_stats <- summary(model)
  c(
    gene = gene,
    HR = summary_stats$coef[1, "exp(coef)"],
    pvalue = summary_stats$coef[1, "Pr(>|z|)"]
  )
})

# Convert to data frame, get q-values and tidy up
cox_results <- do.call(rbind, results)
cox_results <- as.data.frame(cox_results, stringsAsFactors = FALSE)
cox_results$gene <- as.numeric(cox_results$gene)
cox_results$HR <- as.numeric(cox_results$HR)
cox_results$pvalue <- as.numeric(cox_results$pvalue)
qobj <- qvalue(cox_results$pvalue)
cox_results$qvalue <- qobj$qvalues
rownames(cox_results) <- NULL

# Create a dataframe to map Entrez ID to Hugo symbols
entrez_to_hugo <- rna_data[, c("Entrez_Gene_Id", "Hugo_Symbol")]
entrez_to_hugo <- entrez_to_hugo[!duplicated(entrez_to_hugo$Entrez_Gene_Id), ]

# Join the Cox results and the gene mapping dataframe
cox_results <- merge(cox_results, entrez_to_hugo, by.x = "gene", by.y = "Entrez_Gene_Id", all.x = TRUE)

# Rename
cox_results <- cox_results[, c("Hugo_Symbol", "gene", "HR", "pvalue", "qvalue")]
names(cox_results)[2] <- "Entrez_Gene_Id"

cox_results_ordered <- cox_results[order(cox_results$qvalue, decreasing = FALSE), ]

# filtered_cox_results <- cox_results[cox_results$qvalue < 0.2, ]
# cox_results_ordered <- filtered_cox_results[order(filtered_cox_results$HR, decreasing = TRUE), ]

write.csv(cox_results_ordered, "gbm_tcga_pub2013_cox_results.csv", row.names = FALSE)




