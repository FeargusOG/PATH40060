ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

print_banner <- function(section_name) {
  cat("\n")
  cat(strrep("=", 50), "\n")
  cat("SECTION: ", section_name, "\n")
  cat(strrep("=", 50), "\n\n")
}

get_worker_count <- function() {
  ## Use 1 less cores than the total, but ensure at least 1 worker
  return(max(1, parallel::detectCores() - 1))
}

# Define amplification threshold
amplification_threshold <- 0

# Define the standard significnce threshold
SIGNIFICANCE_LEVEL <- 0.05

# ERBB2 Amplification Status Factor Constants
ERBB2_Factor_Levels <- c(0, 1)
ERBB2_Factor_Labels <- c("ERBB2-", "ERBB2+")
ERBB2_Negative_Colour <- "darkblue"
ERBB2_Positive_Colour <- "red"