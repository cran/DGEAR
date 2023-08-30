#' Function for Half's-T-Test Analysis
#'
#' @param con A data frame or matrix containing the expression values for the control.
#' @param exp A data frame or matrix containing the expression values for the experiment.
#' @param alpha Value of significance level ranging from 0 to 1 (default = 0.05 states 5 \% significance).
#' @param FC An array  or list containing fold change values for each gene, calculated by
#'
#' @return A data frame containing values for statistic score, p-values etc for each gene being tested.
#'
#' @export
#'
#' @importFrom stats sd
#' @importFrom stats pt
#' @examples
#' library(DGEAR)
#' data("gene_exp_data")
#' data = read_and_preprocess_data(datafile = gene_exp_data, con1=1,con2=10,exp1=11,exp2=20)
#' perform_h_test(con= data$con, exp= data$exp, FC = data$FC)
perform_h_test <- function(con, exp, alpha= 0.05, FC){
  s0 = apply(con, 1, sd)
  sample_size = sqrt((1/ncol(exp))+(1/ncol(con)))
  half_t_stat = data.frame(rownames(con), (FC)/(s0*sample_size))
  colnames(half_t_stat)= c("ID","statistic")
  for (i in 1:nrow(half_t_stat)) {
    half_t_stat$p.value = 2*pt(abs(half_t_stat$statistic), df=nrow(con)-1, lower.tail = F)
  }
  half_t_stat$fdr[half_t_stat$p.value<=(alpha/nrow(half_t_stat))*seq(length=nrow(half_t_stat))] <- 1
  return(half_t_stat)
}
