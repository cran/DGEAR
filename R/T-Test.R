#' Function for t-Test Analysis
#'
#' @param con A data frame or matrix containing the expression values for the control.
#' @param exp A data frame or matrix containing the expression values for the experiment.
#' @param alpha Value of significance level ranging from 0 to 1 (default = 0.05 states 5 \% significance).
#'
#' @return A data frame containing values for statistic score, p-values etc for each gene being tested.
#'
#' @export
#' @importFrom stats t.test
#' @examples
#' library(DGEAR)
#' data("gene_exp_data")
#' data = read_and_preprocess_data(datafile = gene_exp_data, con1=1,con2=10,exp1=11,exp2=20)
#' perform_t_test(con= data$con, exp= data$exp)
perform_t_test <- function(con, exp, alpha = 0.05) {
  p.value = NULL
  statistic = NULL
  for(i in seq_len(nrow(con))) {
    x = con[i,]
    y = exp[i,]
    t = t.test(x, y)
    p.value[i] = t$p.value
    statistic[i] = t$statistic
  }

  t_stat = cbind.data.frame(ID=row.names(con),p.value,statistic)
  t_stat$fdr[t_stat$p.value<=(alpha/nrow(t_stat))*seq(length=nrow(t_stat))] <- 1
  return(t_stat)
}
