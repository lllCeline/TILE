#' wilcoxon test
#' @description
#' For comparison of two groups,non-parametric two-sided Wilcoxon-rank sum test was used
#' @param score_matrix a data matrix with rows labeling features immune cell type,and columns labeling samples.
#' @param groupinfo a matrix with two columns:'sample' column labeling sampleid,'groups' column labeling group information.
#' when @param paired was setted as'TRUE',@param groupinfo need to provide three columns,the added column 'patient' labeling the paired relation of samples
#' @param paired a logical indicating whether you want a paired wilcoxon test, is a non-parametric alternative to paired t-test used to compare paired data. It's used when your data are not normally distributed,default is "FALSE"
#' @import rstatix
#' @examples
#' result <- wilcox(score_matrix = cell_enrich,groupinfo = demo_groupinfo,paired = F)
#' #paired sample
#' result <- wilcox(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
#' @export

wilcox <- function(score_matrix,groupinfo,paired = FALSE){
  input <- data.frame(score_matrix,check.names = FALSE)
  condition <- data.frame(groupinfo,check.names = FALSE)
  compare <- unique(condition$groups)
  if (TRUE %in% (condition$sample %in% colnames(input))){
    input <- data.frame(t(input),check.names = FALSE)
  }
  if (length(compare) == 2){
    group1 <- compare[1]
    group2 <- compare[2]
    condition <- condition %>%
      group_by(groups) %>%
      filter(n_distinct(sample) > 1) %>%
      ungroup() %>%
      as.data.frame()
    input$sample <- rownames(input)
    input <- input[which(input$sample %in% condition$sample),]
    input <- reshape2::melt(input,variable.name = "cell_type")
    input$group <- condition[match(input$sample,condition$sample,nomatch = 0),2]
    if (paired){
      input$patient <- condition[match(input$sample,condition$sample,nomatch = 0),3]
      input <- input[order(input$group,input$patient),]
      stat.test <-  input %>%
        group_by(cell_type) %>%
        wilcox_test(value ~ group,paired = TRUE) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj")
    }else{
      stat.test <-  input %>%
        group_by(cell_type) %>%
        wilcox_test(value ~ group) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj")
    }
    return(stat.test)
  }else{
    print(compare)
    stop("There are more than two groups")
  }
}
