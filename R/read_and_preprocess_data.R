#' Function to read data and perform initial pre-processing
#'
#' @param datafile A matrix or data frame containing gene expression data
#' @param con1 Starting column of the control of the expression data
#' @param con2 Ending column of the control of the expression data
#' @param exp1 Starting column of the experiment of the expression data
#' @param exp2 Ending column of the experiment of the expression data
#' @param alpha Value of significance level ranging from 0 to 1 (0.05 states 5 \% significance)(Default = 0.05).
#' @param votting_cutoff A numeric value serves as Majority voting (Default = 2)
#'
#' @return A large list containing the data file and the input values
#' @export
#'
#' @examples
#' data("gene_exp_data")
#' read_and_preprocess_data(datafile = gene_exp_data, con1=1,con2=10,exp1=11,exp2=20)
read_and_preprocess_data <- function(datafile, con1, con2, exp1, exp2, alpha = 0.05, votting_cutoff = 2) {
  datafile <- datafile

  con <- datafile[, con1:con2]
  exp <- datafile[, exp1:exp2]

  FC <- rowMeans(exp)-rowMeans(con)

  alpha <- alpha

  # alpha = as.numeric(
  #   readline(prompt = "Specify Alpha value : "))

  votting_cutoff <-votting_cutoff

  # votting_cutoff = as.numeric(
  #   readline(prompt = "Specify votting cutoff(not morethan 5) : "))

  large_list = list(datafile = datafile, con1 = con1, con2 = con2, exp1 = exp1, exp2 = exp2,
       con = con, exp = exp, FC = FC, alpha = alpha, votting_cutoff = votting_cutoff)
  return(large_list)

}

