# Define helper functions and constants
source("R/00_utilities.r")

# 1. Download the dataset on: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
## Performed outside of this script....

# 2. Untar the folder and extract the files.
## Performed outside of this script....

###################################################################
# Parts 3-7:                                                      #
#   - Read the RNA-seq file: data_mrna_seq_v2_rsem.txt            #
#   - Read the Patient Data file: data_clinical_patient.txt       #
#   - Read the Copy Number Aberrations Data: data_cna.txt         #
#   - Match the RNA-seq patient ids with the CNA ids              #
#     and the Patient Data ids.                                   #
#   - Create metadata using the CNA level of ERBB2+               #
#     (greater than 0 means amplified).                           #
###################################################################
print_banner("Preprocessing...")
source("R/01_preprocessing.r")

data <- read_data("data/gbm_tcga_pub2013/data_clinical_patient.txt",
                  "data/gbm_tcga_pub2013/data_mrna_seq_v2_rsem.txt",
                  "data/gbm_tcga_pub2013/data_cna.txt")

## The patients contained in RNASeq and CNA datasets don't match fully; because
## we want to examine the relationship between patients in these two datasets,
## lets drop (the very few) patients from the datasets that aren't in both.
common_cols <- intersect(colnames(data$rnaseq), colnames(data$cna))
data$rnaseq <- data$rnaseq[, common_cols]
data$cna <- data$cna[, common_cols]

meta <- create_amp_metadata(data$rnaseq, data$cna, "EGFR")

###################################################################
# Parts 8-9:                                                      #
#   - Normalize data using DESeq2                                 #
#   - Obtain Differentially Expressed Genes.                      #
###################################################################
print_banner("DEG Analysis...")
source("R/02_DEG.R")
dds <- run_deg_analysis(meta$assay, meta$data, meta$amp_colname)
dds_res <- results(dds)

# Print the top 10 differentially expressed genes by Fold Change.
# The genes listed will have either a positive or negative value for
# `log2FoldChange` depending on whether they are over- or
# under-expressed in the presence of ERBB2+, respectively.
print(dds_res[order(abs(dds_res$log2FoldChange), decreasing = TRUE)[1:10],])

# Get significant, over- and under-expressed DEGs, as a list of symbols
degs <- split_degs(dds_res)

###################################################################
# Part 10:                                                        #
#   - Perform a Pathway Enrichment Analysis                       #
###################################################################
print_banner("Pathway Enrichment...")
source("R/03_pathway_enrichment.R")
# So from our above DESeq analysis, we have found a bunch of genes that are
# either over or under expressed. Now with enrichment analysis we can determine
# what biological processes are invovled in these genes and how statistically
# significant the results are. For example, ERBB2 amplification might cause a
# bunch of other genes related to the cell cycle to become overexpressed. This
# enrichment analysis would say "we're seeing a lot of cell cycle processes
# being impacted in patients with ERBB2 amplification; more than would be
# statistically likely by chance."
run_pathway_enrichment(degs$overx, "Over Expressed", q_thresh = SIGNIFICANCE_LEVEL)
run_pathway_enrichment(degs$underx, "Under Expressed", q_thresh = SIGNIFICANCE_LEVEL)

###################################################################
# Part 11:                                                        #
#   - Get the variance stabilised transformed expression values.  #
###################################################################
print_banner("Variance Stabilisation...")
vsd <- vst(dds)

###################################################################
# Part 12:                                                        #
#   - With the vst values obtain a PCA plot and a heatmap.        #
###################################################################
print_banner("PCA and HeatMap...")
source("R/04_pca_heatmap.R")

###################################################################
# Part 13:                                                        #
#   - With the vst values of the DE genes generate an             #
#     overall survival model using the glmnet package.            #
###################################################################
print_banner("Cox Regression Survival Model...")
source("R/05_cox_regress.R")

X_Y <- prepare_cox_reg_x_y(vsd, data$patient, degs$significant)

X_Y_split <- split_training_data(X_Y$X, X_Y$Y)

cox_models <- fit_cox_models(X_Y_split$X_train, X_Y_split$Y_train)

cox_pred <- make_cox_predictions(cox_models$cox_cv_model, X_Y_split$X_test)

surv_curve <- fit_survival_curves(X_Y_split$Y_test[, "time"],
                                  X_Y_split$Y_test[, "status"], 
                                  ifelse(cox_pred > median(cox_pred),
                                         "High Risk", "Low Risk"))
# Plot the survival curve
ggsurvplot(
  surv_curve$surv, 
  data = surv_curve$data, 
  conf.int = TRUE, 
  title = "Survival Analysis by ERBB2+ Differential Gene Expression", 
  xlab = "Time (Months)"
)


###################################################################
#           Pathway Enrichment for our Cox Model Genes            #
###################################################################
print_banner("Pathway Enrichment for our Cox Model Genes...")

# Extract the genes our model relied on
cox_model_genes <- extract_cox_model_genes(cox_models$cox_cv_model)
# Extract the top 10 genes for high risk (positive coefficients)
top_high_risk_genes <- head(cox_model_genes, 10)
# Extract the top 10 genes for low risk (negative coefficients)
top_low_risk_genes <- tail(cox_model_genes, 10)
print(top_high_risk_genes)
print(top_low_risk_genes)

# Split into risk and protective groups
gene_coefficients <- as.vector(cox_model_genes)
names(gene_coefficients) <- rownames(cox_model_genes)
risk_genes <- names(gene_coefficients[gene_coefficients > 0])
protective_genes <- names(gene_coefficients[gene_coefficients < 0])

# Check the counts
cat("Number of Risk Genes: ", length(risk_genes), "\n")
cat("Number of Protective Genes: ", length(protective_genes), "\n")

run_pathway_enrichment(risk_genes, "Cox Model - Risk Genes")
run_pathway_enrichment(protective_genes, "Cox Model - Protective Genes")


###################################################################
#                     Xu et al. Investigations                    #
###################################################################
print_banner("Xu et al. Investigations...")

# Define the list of genes to check
xu_genes <- c("ACTR6", "C2orf76", "DIO2", "DCXR", "NDUFA8", "SULT1A2", "AQP3")

# Subset the matrix to include only rows matching the xu_genes
matching_genes <- xu_genes %in% rownames(cox_model_genes)

# Print the results
print(matching_genes)

# I din't find any of their genes being used by my model... lets runs the
# survival analysis again with their genes to see if I find them predictive.
Xu_X_Y <- prepare_cox_reg_x_y(vsd, data$patient, xu_genes)
Xu_X_Y_split <- split_training_data(Xu_X_Y$X, Xu_X_Y$Y)

# The model is no good! Basically 0.5 i.e. random chance...
Xu_cox_models <- fit_cox_models(Xu_X_Y_split$X_train, Xu_X_Y_split$Y_train)


###################################################################
#                       FGFR Investigation                        #
###################################################################
print_banner("FGFR Investigation...")

# Process CNA amplifications and generate contingency table and Fisher's test
process_cna_comparison <- function(data_cna, reference_gene, target_gene, amp_thresh) {
  # Process the CNA data for both the reference and target genes
  cna_reference <- extract_cna_gene(data_cna, reference_gene)
  cna_target <- extract_cna_gene(data_cna, target_gene)
  
  # Identify reference gene positive and negative patients
  reference_positive <- colnames(cna_reference)[cna_reference[1, ] > amp_thresh]
  reference_negative <- colnames(cna_reference)[cna_reference[1, ] <= amp_thresh]
  
  # Check amplifications of the target gene in reference gene positive and negative patients
  target_amplified_in_reference_pos <- cna_target[1, reference_positive] > amp_thresh
  target_amplified_in_reference_neg <- cna_target[1, reference_negative] > amp_thresh
  
  # Create contingency table
  contingency_table <- table(
    Reference = c(rep("Positive", length(reference_positive)), rep("Negative", length(reference_negative))),
    Target_Amplified = c(target_amplified_in_reference_pos, target_amplified_in_reference_neg)
  )
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  # Print the results
  cat("\nResults for comparison of", target_gene, "against", reference_gene, "\n\n")
  print(contingency_table)
  print(fisher_test)
}

process_cna_comparison(data$cna, "ERBB2", "FGFR2", amplification_threshold)
process_cna_comparison(data$cna, "ERBB2", "FGFR3", amplification_threshold)
process_cna_comparison(data$cna, "ERBB2", "FGFR4", amplification_threshold)

