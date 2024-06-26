#' Differential Gene Expression Analysis with R
#'
#' @description Main function which incorporates results from five statistical models and detects DEGs through majority voting.
#'
#' @param con1 Starting column of the control of the expression data
#' @param con2 Ending column of the control of the expression data
#' @param exp1 Starting column of the experiment of the expression data
#' @param exp2 Ending column of the experiment of the expression data
#' @param dataframe A matrix containing the gene expression data
#' @param alpha Value of significance level ranging from 0 to 1 (0.05 states 5 \% significance)(Default = 0.05).
#' @param votting_cutoff A numeric value serves as Majority voting (Default = 2)
#'
#' @return A matrix containing Differentially Expressed Genes(DEGs) detected
#'
#' @export
#' @details
#' To use this tool the necessary parameters are con1 = Control start column, con2 = Control end column, exp1 = Experiment start column, exp2 = Experiment end column, alpha = Value of significance level, voting_cutoff = Majority voting value (not more than 5, since there are 5 statistical methods which take part in the majority voting)
#' @importFrom utils read.delim
#' @importFrom stats oneway.test
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @importFrom stats sd
#' @importFrom stats pt
#' @importFrom DescTools DunnettTest
#' @importFrom stats p.adjust
#' @examples
#' library(DGEAR)
#' data("gene_exp_data")
#' DGEAR(dataframe = gene_exp_data, con1 = 1, con2 = 10,
#'   exp1 = 11, exp2 = 20, alpha = 0.05, votting_cutoff = 2)

DGEAR <- function(dataframe, con1, con2, exp1, exp2, alpha, votting_cutoff){
  con= dataframe[,eval(con1:con2)]
  exp= dataframe[,eval(exp1:exp2)]
  con_m = rowMeans(con)
  exp_m = rowMeans(exp)

  #######____t-test____
  p.value = NULL
  statistic = NULL
  for(i in 1 : nrow(dataframe)) {
    x = con[i,]
    y = exp[i,]
    t = t.test(x, y)
    p.value[i] = t$p.value
    statistic[i] = t$statistic
  }

  t_stat = cbind.data.frame(ID=row.names(dataframe),p.value,statistic)
  #t_stat$fdr[t_stat$p.value<=(alpha/nrow(t_stat))*seq(length=nrow(t_stat))] <- 1
  t_stat$BH = p.adjust(t_stat$p.value, method = "BH")
  t_stat$fdr[t_stat$BH<=alpha]<- 1

  #######____anova-oneway-test____
  p.value = NULL
  statistic = NULL
  group = factor(c(rep("con",length(eval(con1:con2))),rep("tre",length(eval(exp1:exp2)))))
  dataframe=t(dataframe)
  for(i in 1 : ncol(dataframe)) {
    x = dataframe[,i]
    o = oneway.test(x~group)
    p.value[i] = o$p.value
    statistic[i] = o$statistic
  }
  dataframe=t(dataframe)
  o_stat = cbind.data.frame(ID=row.names(dataframe),p.value,statistic)
  #o_stat$fdr[o_stat$p.value<=(alpha/nrow(o_stat))*seq(length=nrow(o_stat))] <- 1
  o_stat$BH = p.adjust(o_stat$p.value, method = "BH")
  o_stat$fdr[o_stat$BH<=alpha]<- 1

  #######____Dunnett's test____
  # library(DescTools)
  requireNamespace("DescTools")
  p.value = NULL
  statistic = NULL
  group = factor(c(rep("con",length(eval(con1:con2))),rep("tre",length(eval(exp1:exp2)))))
  dataframe=t(dataframe)
  for(i in 1 : ncol(dataframe)) {
    x = dataframe[,i]
    d = DunnettTest(x,group)
    p.value[i] = d$con[,4]

  }
  dataframe=t(dataframe)
  d_stat = cbind.data.frame(ID=row.names(dataframe),p.value)
  #d_stat$fdr[d_stat$p.value<=(alpha/nrow(d_stat))*seq(length=nrow(d_stat))] <- 1
  d_stat$BH = p.adjust(d_stat$p.value, method = "BH")
  d_stat$fdr[d_stat$BH<=alpha]<- 1


  #######____Half's t-test____
  s0 = apply(con, 1, sd)
  sample_size = sqrt((1/ncol(exp))+(1/ncol(con)))
  half_t_stat = data.frame((exp_m-con_m)/(s0*sample_size))
  colnames(half_t_stat)= "statistic"
  for (i in 1:nrow(half_t_stat)) {
    half_t_stat$p.value = 2*pt(abs(half_t_stat$statistic), df=nrow(con)-1, lower.tail = F)
  }
  #half_t_stat$fdr[half_t_stat$p.value<=(alpha/nrow(half_t_stat))*seq(length=nrow(half_t_stat))] <- 1
  half_t_stat$BH = p.adjust(half_t_stat$p.value, method = "BH")
  half_t_stat$fdr[half_t_stat$BH<=alpha]<- 1



  #######____Wilcox/Mann-wheitneyU-test____
  p.value = NULL
  statistic = NULL
  for(i in 1 : nrow(dataframe)) {
    x = as.numeric(con[i,])
    y = as.numeric(exp[i,])
    u = wilcox.test(x, y, paired = F)
    p.value[i] = u$p.value
    statistic[i] = u$statistic
  }

  u_stat = cbind.data.frame(ID=row.names(dataframe),p.value,statistic)
  #u_stat$fdr[u_stat$p.value<=(alpha/nrow(u_stat))*seq(length=nrow(u_stat))] <- 1
  u_stat$BH = p.adjust(u_stat$p.value, method = "BH")
  u_stat$fdr[u_stat$BH<=alpha]<- 1



  ######### ____ENSEMBL_____
  df = cbind.data.frame(as.numeric(t_stat$fdr),
                        as.numeric(o_stat$fdr),
                        as.numeric(d_stat$fdr),
                        as.numeric(half_t_stat$fdr),
                        as.numeric(u_stat$fdr))
  row.names(df)= t_stat$ID
  df[is.na(df)] <- 0


  DF = data.frame(DEGs = row.names(df[rowSums(df) >= votting_cutoff,]))

  return(DF)



}
