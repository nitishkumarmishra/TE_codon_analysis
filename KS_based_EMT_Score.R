# By Nitis 11/01/2024
# This code is based on these two papers-
# https://www.embopress.org/doi/full/10.15252/emmm.201404208
# https://www.pnas.org/doi/10.1073/pnas.2102050118
# This script calculates EMT (Epithelial-Mesenchymal Transition) scores using the Kolmogorov-Smirnov (K-S) test.
# The EMT score is determined by comparing the cumulative distribution functions (CDFs) of epithelial and mesenchymal gene signatures.
# The score ranges from -1 to +1, where a positive score indicates a mesenchymal phenotype and a negative score indicates an epithelial phenotype.

#I can use Epithelial and mesenchymal marker genes, which I have from literature, rather than based on top epithelial genes and then calculate correlation.
#Select genes which have high positive and negative correlation and expand epithelial and mesenchymal gene list.


# Load necessary library
library(stats)

# Example data frames (replace with your actual data)
epithelial_data <- data.frame(
  sample1 = rnorm(100, mean = 5, sd = 2),
  sample2 = rnorm(100, mean = 6, sd = 2)
)

mesenchymal_data <- data.frame(
  sample1 = rnorm(100, mean = 7, sd = 2),
  sample2 = rnorm(100, mean = 8, sd = 2)
)

# Function to calculate EMT score using K-S test
calculate_emt_score <- function(epithelial_data, mesenchymal_data) {
  emt_scores <- sapply(1:ncol(epithelial_data), function(i) {
    epithelial_sample <- epithelial_data[, i]
    mesenchymal_sample <- mesenchymal_data[, i]
    
    # Calculate the K-S test statistic
    ks_test <- ks.test(epithelial_sample, mesenchymal_sample)
    ks_statistic <- ks_test$statistic
    
    # Determine the direction of the score
    if (mean(mesenchymal_sample) > mean(epithelial_sample)) {
      return(ks_statistic)
    } else {
      return(-ks_statistic)
    }
  })
  return(emt_scores)
}

# Calculate EMT scores
emt_scores <- calculate_emt_score(epithelial_data, mesenchymal_data)

# Print EMT scores
print(emt_scores)
