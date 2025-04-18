library(dplyr)

read_data <- function(patient, rnaseq, cna) {
  ## Skip first 4 rows as they don't contain observations.
  data_patient <- read.delim(patient, skip = 4)
  data_rnaseq <- read.delim(rnaseq)
  data_cna <- read.delim(cna)
  
  return(list(patient = data_patient,
              rnaseq = data_rnaseq,
              cna = data_cna))
}

extract_cna_gene <- function(data_cna, gene_symbol) {
  gene_cna <- data_cna[data_cna$Hugo_Symbol == gene_symbol, ] |>
    dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id)
  rownames(gene_cna) <- gene_symbol
  return(gene_cna)
}

create_amp_metadata <- function(data_rnaseq, data_cna, gene, amp_thresh = 0.05) {
  cna_extracted_gene <- extract_cna_gene(data_cna, gene)

  ## Extract assay data, skipping the first two columns (gene identifiers).
  assay <- round(as.matrix(data_rnaseq[, -c(1, 2)]))
  rownames(assay) <- data_rnaseq[, 1]

  ## Create an empty matrix to store our metadata.
  metadata <- matrix(0, ncol = 1, nrow = ncol(assay))

  ## Loop over each column in the assay dataset.
  for (i in seq_len(ncol(assay))) {
    ## Find the patient with the matching ID in the CNA dataset.
    patient_col <- cna_extracted_gene |>
      dplyr::select(colnames(assay)[i])
    ## Check the amplification level
    metadata[i, 1] <- ifelse(as.numeric(patient_col[1,1]) > amp_thresh, 1, 0)
  }

  ## Replace NA values in metadata with 0
  metadata[is.na(metadata)] <- 0
  amp_colname <- paste0(gene,"_Amp")
  colnames(metadata) <- c(amp_colname)
  return(list(assay = assay,
              data = metadata,
              amp_colname = amp_colname))
}
