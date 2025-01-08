#' Title
#'
#' @param data A dataframe where rows represent unique genes and columns
#' represent samples. Gene identifiers (rownames) can be provided as HUGO gene
#' symbols.RNA-seq data should be log-transformed,
#' for example using log2(normalized counts + 1).
#' @param cor_cut  A numeric value specifying the minimum Pearson correlation
#' threshold for classification.Samples with a maximum correlation below this
#' threshold will remain unclassified, and their classification results will be
#' set to NA. Default value is 0.1.
#'
#' @return A dataframe containing classification results for each sample.
#' The `consensus_cluster` column provides the predicted consensus class label
#' for each sample. Additional columns include the Pearson correlation values
#' for each sample with each centroid, the `normalized_difference`
#' (a measure of how distinct the sample is from other classes),
#' and the `p_value` associated with the correlation to the nearest centroid.
#' Samples with a maximum correlation below the `cor_cut` threshold will have
#' NA in the `consensus_cluster` column.
#' @export
#'
#' @examples
#' # Example usage:
#' # `data` is a gene expression matrix:
#' data(tcgadat)
#' classify_nmibc(dejongA, cor_cut = 0.1)
classify_nmibc <- function(data, cor_cut = 0.1) {
  # Load the centroids matrix from the NMIBCClassifier package
  centroids <- NMIBCClassifier::centroids
 
  # Identify common genes between the input data and the centroids
  common_genes <- intersect(rownames(centroids), rownames(data))
 
  # Check if the number of common genes is less than half of the centroid genes
  if (length(common_genes) < nrow(centroids) / 2) {
    warning("The number of common genes is less than half of the centroid genes. The results may be unstable.")
  }
  
  # Define a helper function to calculate Pearson correlation
  pearson_correlation <- function(x, y) {
    cor(x, y, method = "pearson")
  }

  results <- t(apply(data[common_genes, , drop = FALSE], 2, function(sample) {

    # Compute Pearson correlations between the sample and each centroid
    correlations <- apply(centroids[common_genes, , drop = FALSE],
                          2,
                          pearson_correlation,
                          y = sample)

    # Identify the class with the highest correlation
    max_correlation <- max(correlations)
    classification <- colnames(centroids)[which.max(correlations)]
    # Ensure the class is a column name of centroids

    # Identify the second highest correlation
    second_max_correlation <- sort(correlations, decreasing = TRUE)[2]

    # Calculate the normalized difference
    normalized_difference <-
      (max_correlation - second_max_correlation) / max_correlation

    # Compute the p-value for the correlation with the nearest centroid
    nearest_centroid <- centroids[common_genes, classification]
    # Ensure correct indexing
    cor_pval <- cor.test(sample, nearest_centroid)$p.value

    # If the maximum correlation < threshold, set classification to NA
    if (max_correlation < cor_cut) {
      classification <- NA
    }

    # Return classification results, correlations, normalized difference and p
    c(
      classification = classification,
      correlations,
      normalized_difference = normalized_difference,
      p_value = cor_pval
    )
  }))

  # Convert the results to a dataframe
  results_df <- as.data.frame(results)

  colnames(results_df) <- c(
    "consensus_cluster",
    paste0("correlation_to_consensus_cluster_", seq_len(ncol(centroids))),
    "normalized_difference",
    "p_value"
  )

  return(results_df)
}
