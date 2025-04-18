## Ensure BiocManager and DESeq2 are installed
ensure_package("BiocManager")
ensure_package("BiocParallel") # DESeq2 was slow! Let's use more cores....
ensure_package("DESeq2")
library(BiocParallel)
library(DESeq2)

run_deg_analysis <- function(assay, metadata, gene_amp) {
  # 8 - Normalize data using DESeq2
  ## Prepare RNASeq data for DESeq2
  assay[is.na(assay)] <- 0  # Replace NA values with zeros
  assay[assay < 0] <- 0     # Set negative values to zero
  
  ## Filter out genes with low expression or expressed in few samples; this
  ## protects our analysis from spurious results due to noise etc.
  smallest_group_size <- 3
  keep <- rowSums(assay >= 10) >= smallest_group_size
  assay <- assay[keep, ]
  
  # Use as.formula() to dynamically create the formula
  design_formula <- as.formula(paste("~", gene_amp))
  
  ## Create the differential expression analysis dataset.
  ##   assay:    each row is a gene and each column a patient.
  ##   metadata: each row is a patient and the single column is a grouping factor.
  ##
  ## Therefore, the rows in metadata must match the columns in countData in both 
  ## number and order!
  ##
  ## Basically, what this function does is say "here is a bunch of patients with 
  ## their genome data. Here is a grouping factor (could be age, could be ERB22+ 
  ## expression, whatever). What patterns can be see in their genome exoression 
  ## data from this grouping factor?"
  dds <- DESeqDataSetFromMatrix(countData = assay,
                                colData = metadata,
                                design = design_formula)
  
  ## Run DESeq2 differential expression analysis
  register(MulticoreParam(workers = get_worker_count()))
  dds <- DESeq(dds, parallel = TRUE)
  return(dds)
}

split_degs <- function(dds_res, p_thresh = 0.05) {
  degs_significant = dds_res[dds_res$padj < p_thresh,]
  degs_overx = rownames(degs_significant[degs_significant$log2FoldChange>0,])
  degs_underx = rownames(degs_significant[degs_significant$log2FoldChange<0,])
  
  # Print some stats
  cat("Number of significant DEGs found:", nrow(degs_significant), 
      "out of a total of", nrow(dds_res), "genes (",
      round(100 * nrow(degs_significant) / nrow(dds_res), 2), "%)\n")
  cat("  - Over-expressed genes:", length(degs_overx), "\n")
  cat("  - Under-expressed genes:", length(degs_underx), "\n")
  
  return(list(significant = rownames(degs_significant),
              overx = degs_overx,
              underx = degs_underx))
}
