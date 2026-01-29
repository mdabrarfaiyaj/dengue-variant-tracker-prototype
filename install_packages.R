#!/usr/bin/env Rscript
# Installation script for Dengue Variant Tracker Dashboard
# Run this once to install all required R packages

cat("==========================================\n")
cat("Installing R Packages for Dengue Tracker\n")
cat("==========================================\n\n")

# Function to check and install packages
install_if_missing <- function(pkg, from = "CRAN") {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s from %s...\n", pkg, from))
    if (from == "CRAN") {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    } else if (from == "Bioconductor") {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    cat(sprintf("✓ %s installed successfully\n\n", pkg))
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# Install CRAN packages
cat("\n[1/2] Installing CRAN packages...\n")
cat("------------------------------------\n")
cran_packages <- c(
  "shiny",
  "shinydashboard", 
  "ggplot2",
  "dplyr",
  "tidyr",
  "plotly",
  "DT"
)

for (pkg in cran_packages) {
  install_if_missing(pkg, "CRAN")
}

# Install Bioconductor packages
cat("\n[2/2] Installing Bioconductor packages...\n")
cat("------------------------------------\n")
bioc_packages <- c(
  "Biostrings",
  "ShortRead",
  "BSgenome"
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, "Bioconductor")
}

# Verify installation
cat("\n==========================================\n")
cat("Verifying Installation...\n")
cat("==========================================\n")

all_packages <- c(cran_packages, bioc_packages)
missing <- character(0)

for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing <- c(missing, pkg)
    cat(sprintf("✗ %s - NOT FOUND\n", pkg))
  } else {
    cat(sprintf("✓ %s - OK\n", pkg))
  }
}

cat("\n==========================================\n")
if (length(missing) == 0) {
  cat("SUCCESS! All packages installed correctly.\n")
  cat("\nNext steps:\n")
  cat("  1. Run: ./download_data.sh\n")
  cat("  2. Run: Rscript qc_analysis.R\n")
  cat("  3. Run: Rscript -e \"shiny::runApp('app.R')\"\n")
} else {
  cat("WARNING: Some packages failed to install:\n")
  cat(paste("  -", missing, collapse = "\n"))
  cat("\n\nPlease install them manually:\n")
  cat("  install.packages(c('", paste(missing, collapse = "', '"), "'))\n", sep = "")
}
cat("==========================================\n")
