# Set working dir
setwd("working_dir_path")


tmpdir <- "tmp_dir_path"
options(tempdir = tmpdir)

# Check libraries
required_packages <- c("LDlinkR")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "https://cloud.r-project.org", lib = "~/R/library")
  }
}

# Load lib
library(LDlinkR, lib.loc = "~/R/library")

# Add access token
access_token <- "access_token_here"

# Load Data
df <- read.table("cardio_gwas_SNPs", sep = "\t", header = TRUE)
rsids_df <- df$rsID
rsids_df <- as.data.frame(rsids_df)

# Variable initialization
ld_results <- NULL
error_log <- NULL

# Set splitting size (chunks of 10)
chunk_size <- 10
num_chunks <- ceiling(nrow(df) / chunk_size)

# Loop through all chunks
for (i in 1:num_chunks) {
  # Set the start and end of current chunk
  start_index <- (i - 1) * chunk_size + 1
  end_index <- min(i * chunk_size, nrow(df))

  # Extract the rsIDs from chunk
  rsids <- rsids_df[start_index:end_index, 1]

  # Define tissues and populations
  tissues <- c("Heart_Left_Ventricle")
  populations <- c("AFR", "EUR", "SAS", "EAS", "AMR")

  # Run LDexpress
  tryCatch({
    ld_result <- LDexpress(snps = rsids, tissue = tissues, pop = populations,
                            r2d = "r2", r2d_threshold = 0.8, p_threshold = 0.1,
                            win_size = 500000, genome_build = "grch38", token = access_token)

    # Append ld_result to ld_results
    if (is.null(ld_results)) {
      ld_results <- ld_result
    } else {
      ld_results <- rbind(ld_results, ld_result)
    }
  }, error = function(e) {
    # Log errors
    error_log <- c(error_log, paste("Error:", e$message))
    error_log <- c(error_log, paste("Problematic rsIDs:", paste(rsids, collapse = ", ")))
  })
}

# Export 'ld_results'
write.csv(ld_results, file = "ld_results_left_ventricle.csv", row.names = FALSE)

# Export error_log
if (length(error_log) > 0) {
  writeLines(error_log, con = "error_log.txt")
  cat("Error log saved to 'error_log.txt'.\n")
} else {
  cat("No errors encountered.\n")
}