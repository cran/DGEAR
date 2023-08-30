#' Function for Wilcoxon-Mann-Whitney U-Test
#'
#' @param con A data frame or matrix containing the expression values for the control.
#' @param exp A data frame or matrix containing the expression values for the experiment.
#' @param alpha Value of significance level ranging from 0 to 1 (default = 0.05 states 5 \% significance).
#'
#' @return A data frame containing values for statistic score, p-values etc for each gene being tested.
#'
#' @export
#'
#' @importFrom stats wilcox.test
#'
#' @examples
#' library(DGEAR)
#' data("gene_exp_data")
#' data = read_and_preprocess_data(datafile = gene_exp_data, con1=1,con2=10,exp1=11,exp2=20)
#' perform_wilcox_test(con= data$con, exp= data$exp)
perform_wilcox_test <- function(con, exp, alpha = 0.05) {
  p.value = NULL
  statistic = NULL
  for(i in 1 : nrow(con)) {
    x = as.numeric(con[i,])
    y = as.numeric(exp[i,])
    u = wilcox.test(x, y, paired = F)
    p.value[i] = u$p.value
    statistic[i] = u$statistic
  }

  u_stat = cbind.data.frame(ID=row.names(con),p.value,statistic)
  u_stat$fdr[u_stat$p.value<=(alpha/nrow(u_stat))*seq(length=nrow(u_stat))] <- 1
  return(u_stat)
}
