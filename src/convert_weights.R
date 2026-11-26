# Consolidate all cis-expression prediction weights into TSV
# The original approach was going to be to write 3 tsv files for each rds file
# However, the weights_cis_pval_elnet contains 7182 items, leading to 21546 files
# So I opted to consolidate them all with this script instead

# Define directories

input_dir <- "data/weights_cis_pval_elnet"
output_dir <- "data/expression_weights"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Identify all RDS files

rds_files <- list.files(
  path = input_dir,
  pattern = "\\.weights\\.rds$",
  full.names = TRUE
)

total_files <- length(rds_files)

cat("\nGCSC EXPRESSION WEIGHT CONVERSION\n")
cat("RDS files found:", total_files, "\n\n")

if (total_files == 0) {
  stop("No .weights.rds files found.")
}

# Initialize storage and counters

all_weights <- vector("list", total_files)

successful <- 0
failed <- 0

expected <- c(
  "weights_a",
  "weights_b",
  "weights_co"
)

# Process each RDS file

for (i in seq_along(rds_files)) {

  input_file <- rds_files[i]

  file_name <- basename(input_file)

  gene_pair <- sub(
    "\\.weights\\.rds$",
    "",
    file_name
  )

  tryCatch({

    # Read original RDS object

    weights <- readRDS(input_file)

    if (!all(expected %in% names(weights))) {
      stop("Missing expected components.")
    }

    # Extract components

    file_data <- lapply(expected, function(component) {

      x <- weights[[component]]

      if (is.null(names(x)) || length(x) == 0) {
        stop(paste("Invalid SNP identifiers in", component))
      }

      if (anyNA(names(x)) || any(names(x) == "")) {
        stop(paste("Missing SNP identifiers in", component))
      }

      if (!is.numeric(x)) {
        stop(paste("Non-numeric coefficients in", component))
      }

      if (any(!is.finite(x))) {
        stop(paste("Non-finite coefficients in", component))
      }

      data.frame(
        SNP = names(x),
        weight = as.numeric(x),
        component = component,
        gene_pair = gene_pair,
        stringsAsFactors = FALSE
      )
    })

    # Combine the three components for this file

    all_weights[[i]] <- do.call(rbind, file_data)

    successful <- successful + 1

  }, error = function(e) {

    failed <<- failed + 1

    cat(
      "\nFailed:",
      file_name,
      "|",
      conditionMessage(e),
      "\n"
    )
  })

  # Progress update every 500 files

  if (i %% 500 == 0 || i == total_files) {

    cat(
      "Progress:",
      i, "/", total_files,
      "| Successful:", successful,
      "| Failed:", failed,
      "\n"
    )
  }
}

# Combine all successful files

all_weights <- Filter(
  Negate(is.null),
  all_weights
)

if (length(all_weights) == 0) {
  stop("No weight files were successfully processed.")
}

final_weights <- do.call(
  rbind,
  all_weights
)

# Export consolidated TSV

output_file <- file.path(
  output_dir,
  "expression_weights.tsv"
)

write.table(
  final_weights,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Final summary

cat("CONVERSION COMPLETE\n")

cat("Total RDS files:", total_files, "\n")
cat("Successful files:", successful, "\n")
cat("Failed files:", failed, "\n")

cat("Total exported rows:", nrow(final_weights), "\n")

cat(
  "Unique gene pairs:",
  length(unique(final_weights$gene_pair)),
  "\n"
)

cat(
  "Unique SNP identifiers:",
  length(unique(final_weights$SNP)),
  "\n"
)

cat("Output:", output_file, "\n")