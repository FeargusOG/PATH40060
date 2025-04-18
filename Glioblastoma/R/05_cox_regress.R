library(glmnet)
library(survival)
library(survminer)
library(doParallel)

standardise_ids <- function(matrix) {
  # Extract row names
  row_names <- rownames(matrix)
  
  # Standardise IDs
  row_names <- sub("\\.\\d+$", "", row_names) # Remove trailing ".01"
  row_names <- gsub("\\.", "-", row_names)   # Replace '.' with '-'
  
  # Update row names in the matrix
  rownames(matrix) <- row_names
  
  return(matrix)
}

prepare_vsd_matrix <- function(vsd, gene_list) {
  # Subset the vst data for only specified genes
  vsd_genes_matrix <- assay(vsd)[gene_list, ]
  
  # Transpose the matrix to have patients as rows and genes as columns
  vsd_genes_matrix_t <- t(vsd_genes_matrix)
  
  # Standardise the IDs
  vsd_genes_matrix_t <- standardise_ids(vsd_genes_matrix_t)
  
  return(vsd_genes_matrix_t)
}

prepare_patient_survivability <- function(patient_data) {
  ## Set the PATIENT_ID column as the rownames
  Y_patient <- data.frame(
    patient_id = patient_data$PATIENT_ID,
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

prepare_cox_reg_x_y <- function(vsd, patient_data, gene_list) {
  ## Preprocessing for X
  X_vsd <- prepare_vsd_matrix(vsd, gene_list)
  
  ## Preprocessing for Y
  Y_patient <- prepare_patient_survivability(patient_data)
  
  ##
  ## Match X and Y
  ##
  ## The patients in each dataset don't match perfectly;
  ## Find the common rownames between X_vsd and Y_patient
  common_rownames <- intersect(rownames(X_vsd), rownames(Y_patient))
  
  ## Subset X_vsd and Y_patient to only include rows with common rownames
  X <- X_vsd[common_rownames, , drop = FALSE]
  Y_patient <- Y_patient[common_rownames, , drop = FALSE]
  
  # Order both X and Y_patient explicitly by their row names
  X <- X[order(rownames(X)), , drop = FALSE]
  Y_patient <- Y_patient[order(rownames(Y_patient)), , drop = FALSE]
  
  # Ensure both are aligned by their row names after ordering
  stopifnot(all(rownames(X) == rownames(Y_patient)))
  
  ## Create the Surv object
  Y <- Surv(time = Y_patient$time, event = Y_patient$status)
  
  return(list(X = X, Y = Y))
}

split_training_data <- function(X, Y, seed = 1234, split = 0.8) {
  ##
  ## Split training and test data
  ##
  # Set seed for reproducibility
  set.seed(seed)
  
  # Get the indices for the training set (80% of the data)
  train_ind <- sample(1:nrow(X), size = floor(split * nrow(X)))
  
  # Split X into training and test sets
  X_train <- X[train_ind, , drop = FALSE]
  X_test <- X[-train_ind, , drop = FALSE]
  
  # Split Y into training and test sets
  Y_train <- Y[train_ind]
  Y_test <- Y[-train_ind]
  
  # Confirm the dimensions
  cat("Training set: ", nrow(X_train), "samples\n")
  cat("Test set: ", nrow(X_test), "samples\n")
  
  return(list(X_train = X_train,
              Y_train = Y_train,
              X_test = X_test,
              Y_test = Y_test))
}

# SOURCE
# https://glmnet.stanford.edu/articles/Coxnet.html
fit_cox_models <- function(X_train, Y_train, plot_it = TRUE) {
  # Its slow! Let's paralelise it.
  cl <- makeCluster(get_worker_count())
  registerDoParallel(cl)
  
  cox_model <- glmnet(X_train, Y_train, family = "cox", parallel = TRUE)
  cox_cv_model <- cv.glmnet(X_train, Y_train, family = "cox", type.measure = "C", parallel = TRUE)
  
  if(plot_it){
    plot(cox_model)
    plot(cox_cv_model)
  }
  
  # Stop the parelel stuff.
  stopCluster(cl)
  
  return(list(cox_model = cox_model, cox_cv_model = cox_cv_model))
}

make_cox_predictions <- function(cox_cv_model, X_test, lambda = "lambda.min") {
  # Generate risk scores for the test set
  scores <- predict(cox_cv_model, newx = X_test, s = lambda)
  return(as.numeric(scores[, 1]))  # Extract the first (and only) column as a numeric vector
}

# SOURCE
# https://bookdown.org/staedler_n/highdimstats/survival-analysis.html#regularized-cox-regression
fit_survival_curves <- function(Y_time, Y_status, group_by) {
  dat_test <- data.frame(
    time = Y_time, 
    status = Y_status, 
    strat = group_by
  )
  
  # Fit a survival curve using the test data
  return(list(
    surv = survfit(Surv(time, status) ~ strat, data = dat_test), 
    data = dat_test))
}

extract_cox_model_genes <- function(cox_cv_model, lambda = "lambda.min") {
  # Extract coefficients at the optimal lambda as a matrix
  non_zero_genes <- as.matrix(coef(cox_cv_model, s = lambda))
  non_zero_genes <- non_zero_genes[non_zero_genes[, 1] != 0, , drop = FALSE]
  # Sort the non-zero genes by coefficient values
  return(non_zero_genes[order(non_zero_genes[, 1], decreasing = TRUE), , drop = FALSE])
}