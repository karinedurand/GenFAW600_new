#!/usr/bin/env Rscript

# ======================================================================
# Install LEA package (uncomment if not already installed)
# ======================================================================
# install.packages("BiocManager")
# BiocManager::install("LEA", update = TRUE, ask = FALSE)

# Load the LEA package
library(LEA)

# ======================================================================
# Run sNMF analysis (example for K = 2 to 8 — adjust as needed)
# ======================================================================


vcf_file <- "WholeGenome_biallelic_max5missing_pruned.vcf"
##############################################################################################
# Run sNMF with alpha = 100
##############################################################################################

proj_alpha100 <- snmf(vcf_file, 
                      K = 2:8,              # Range of K values to test
                      ploidy = 2,           # Diploid organisms (humans, animals, plants)
                      entropy = TRUE,       # Compute cross-entropy criterion
                      repetitions = 5,     # Number of runs per K for stability
                      iterations = 200,     # Iterations per run (increase for higher accuracy)
                      CPU = 8,              # Number of CPU cores for parallelization
                      project = "new" ,
                      alpha = 100,          # Regularization parameter
                      seed = 123)           # Random seed for reproducibility

# ======================================================================
# Save the sNMF project
# ======================================================================
save(proj_alpha100, file = "WholeGenome_biallelic_max5missing_pruned.RData")

# ======================================================================
# Plot cross-entropy results for alpha = 100
# ======================================================================

# Load the sNMF project (use the same name as above)

project_alpha100 <- load.snmfProject("WholeGenome_biallelic_max5missing_pruned.snmfProject")

# Define the range of K values
Ks <- 2:8

# Compute mean cross-entropy for each K
ce <- sapply(Ks, function(k) mean(cross.entropy(project_alpha100, K = k)))

# Identify the best K (minimum cross-entropy)
best_K<- Ks[which.min(ce)]
cat("Best K (alpha = 100):", best_K, "\n")
# Save the cross-entropy plot
pdf("cross_entropy_alpha100.pdf")
plot(Ks, ce, type = "b", pch = 19, col = "blue",
     xlab = "K (number of clusters)",
     ylab = "Cross-Entropy",
     main = "Cross-Entropy Criterion (sNMF, alpha = 100)")
dev.off()


##############################################################################################
# Plot ancestry barplots with ordered metadata
##############################################################################################

# Load metadata file (edit the filename and separator if needed)
metadata <- read.table("metadata_sNMF.ind", header = TRUE, sep = ",")

# Get Q-matrix for the best K
Q_matrix <- Q(project, K = best_K)

# Add a column for original sample order
metadata$index <- 1:nrow(metadata)

# Ensure that population is treated as a factor
metadata$pop <- factor(metadata$pop)

# Reorder samples by population
order_idx <- order(metadata$pop)
Q_matrix_ordered <- Q_matrix[order_idx, ]
metadata_ordered <- metadata[order_idx, ]

# Save ordered outputs (useful for reproducibility)
write.csv(Q_matrix_ordered, "Q_matrix_ordered_alpha100.csv", row.names = FALSE)
write.csv(metadata_ordered, "metadata_ordered_alpha100.csv", row.names = FALSE)

# Labels (population names)
labels <- metadata_ordered$pop

# Save the ordered admixture barplot for the best K
pdf("admixture_barplot_pop_alpha100_bestK.pdf", width = 12, height = 6)
barplot(t(Q_matrix_ordered),
        col = rainbow(best_K),
        border = NA,
        space = 0,
        main = paste("Genetic Structure - K =", best_K  ),
        names.arg = labels,
        las = 2,               # Rotate labels vertically
        cex.names = 0.6)       # Label size
dev.off()


##############################################################################################
# Plot admixture barplots for all tested K values
##############################################################################################

# Loop over all K values tested in the project
for (k in Ks) {
  Qk <- Q(project, K = k)
  pdf(paste0("admixture_barplot_alpha100_K", k, ".pdf"), width = 12, height = 6)
  barplot(t(Qk[order_idx, ]),
          col = rainbow(k),
          border = NA,
          space = 0,
          main = paste("Genetic Structure - K =", k),
          names.arg = labels,
          las = 2,
          cex.names = 0.6)
  dev.off()
}

