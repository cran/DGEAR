#' Function for Dunnett's Test
#'
#' @param datafile A matrix containing the gene expression data
#' @param con A data frame or matrix containing the expression values for the control.
#' @param exp A data frame or matrix containing the expression values for the experiment.
#' @param alpha Value of significance level ranging from 0 to 1 (default = 0.05 states 5 \% significance).
#'
#' @return A data frame containing values for statistic score, p-values etc for each gene being tested.
#' @export
#' @importFrom DescTools DunnettTest
#'
#' @examples
#' library(DGEAR)
#' data("gene_exp_data")
#' data = read_and_preprocess_data(datafile = gene_exp_data, con1=1,con2=10,exp1=11,exp2=20)
#' perform_dunnett_test(datafile = data$datafile, con= data$con, exp= data$exp)
perform_dunnett_test <- function(datafile, con, exp, alpha =0.05) {
  # library(DescTools)
  requireNamespace("DescTools")
  p.value = NULL
  statistic = NULL
  group = factor(c(rep("con",length(con)),rep("tre",length(exp))))
  dataframe=t(datafile)
  for(i in 1 : ncol(dataframe)) {
    x = dataframe[,i]
    d = DunnettTest(x,group)
    p.value[i] = d$con[,4]

  }
  d_stat = cbind.data.frame(ID=row.names(con),p.value)
  d_stat$fdr[d_stat$p.value<=(alpha/nrow(d_stat))*seq(length=nrow(d_stat))] <- 1
  return(d_stat)
}
