path <- "data/weights_cis_pval_elnet/AARSD1-DDHD2.weights.rds"

weights <- readRDS(path)

for (component in names(weights)) {

  cat("COMPONENT:", component, "\n")

  x <- weights[[component]]

  cat("\nNumber of SNPs:", length(x), "\n")

  cat("\nFirst 10 SNP identifiers:\n")
  print(head(names(x), 10))

  cat("\nFirst 10 coefficients:\n")
  print(head(x, 10))

  cat("\nCoefficient summary:\n")
  print(summary(x))

  cat("\nZero coefficients:", sum(x == 0), "\n")

  cat("Nonzero coefficients:", sum(x != 0), "\n")

  cat("Missing coefficients:", sum(is.na(x)), "\n")
}