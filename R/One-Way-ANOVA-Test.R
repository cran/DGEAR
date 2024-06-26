#' Function for ANOVA One-Way Test
#'
#' @param datafile A matrix containing the gene expression data
#' @param con A data frame or matrix containing the expression values for the control.
#' @param exp A data frame or matrix containing the expression values for the experiment.
#' @param alpha Value of significance level ranging from 0 to 1 (default = 0.05 states 5 \% significance).
#'
#' @return A data frame containing values for statistic score, p-values etc for each gene being tested.
#' @export
#' @importFrom stats p.adjust
#' @importFrom stats na.omit
#' @examples
#' library(DGEAR)
#' data("gene_exp_data")
#' data = read_and_preprocess_data(datafile = gene_exp_data, con1=1,con2=10,exp1=11,exp2=20)
#' perform_anova(datafile = data$datafile, con= data$con, exp= data$exp)
perform_anova <- function(datafile, con, exp, alpha =0.05) {
  p.value = NULL
  statistic = NULL
  group = factor(c(rep("con",length(con)),rep("tre",length(exp))))
  dataframe=t(as.data.frame(datafile))
  for(i in 1 : ncol(dataframe)) {
    x = dataframe[,i]
    o = oneway.test(x~group)
    p.value[i] = o$p.value
    statistic[i] = o$statistic
  }
  o_stat = cbind.data.frame(ID=row.names(datafile),statistic,p.value)
  #o_stat$fdr[o_stat$p.value<=(alpha/nrow(o_stat))*seq(length=nrow(o_stat))] <- 1
  o_stat$BH = p.adjust(o_stat$p.value, method = "BH")
  o_stat$fdr[o_stat$BH<=alpha]<- 1
  
  DEGs=o_stat$ID[o_stat$fdr == 1]
  DEGs= na.omit(DEGs)
  DEGs = as.data.frame(DEGs)
  colnames(DEGs) = "DEGs"
  print(DEGs)
  return(list(Table = o_stat, DEGs = DEGs))
}
