#' A dataset containing gene expression data
#'
#' This dataset contains statistically simulated gene expression data for ease of exercise.
#'
#' @docType data
#' @name gene_exp_data
#' @format A data frame with 10 rows and 20 columns, the columns represents samples, say first 10 columns 1 to 10 being control and 11 to 20 being experiment. Whereas, the rows of the dataset contains genes. First 5 out of 10 genes, gene1-gene5 are the true DEGs as the expression values for the first 10 samples are ~13 times higher than the rest. 
#' @examples
#' # Data will be loaded with lazy loading and can be accessible when needed.
#' data("gene_exp_data")
#' head(gene_exp_data)
"gene_exp_data"
