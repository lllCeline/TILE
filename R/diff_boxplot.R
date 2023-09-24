#' multiple-panels box plot
#' @description
#' specify cell type variables for faceting the plot into multiple panels
#' @param score_matrix a data matrix with rows labeling features immune cell type,and columns labeling samples
#' @param paired a logical indicating whether you want a paired wilcoxon test, which is used to compare the means of two samples when each observation in one sample can be paired with an observation in the other sample,default is "FALSE"
#' @param groupinfo a matrix with two columns:'sample' column labeling sampleid,'groups' column labeling group information.
#' when @param paired was TRUE, @param groupinfo need to provide three columns,the added column 'patient' labeling the paired relation of samples
#' @param fontsize a numeric value change the font size, default is NULL(function can provide the font size automatically)
#' @import ggplot2
#' @import ggpubr
#' @import reshape2
#' @examples
#' p <- boxplot_diff(score_matrix = cell_enrich,groupinfo = demo_groupinfo,paired = F)
#' print(p)
#' #paired sample
#' p <- boxplot_diff(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
#' print(p)
#' @export


boxplot_diff <- function(score_matrix,groupinfo,paired = FALSE,fontsize = NULL){
  data <- data.frame(score_matrix, check.names = FALSE)
  condition <- data.frame(groupinfo,check.names = FALSE)
  compare <- unique(condition$groups)
  if (TRUE %in% (condition$sample %in% colnames(data))){
    data <- data.frame(t(data),check.names = FALSE)
  }
  cell_type <- colnames(data)
  data <- data[which(rownames(data) %in% condition$sample),,drop = FALSE]
  data <- as.data.frame(t(data))
  if (length(compare) == 2){
    group1 <- compare[1]
    group2 <- compare[2]
    data <- subset(data,select = as.character(condition$sample))
    data$Cell <- rownames(data)
    data2 <- melt(data,variable.name = "sample")
    data2$Group <- condition[match(data2$sample,condition$sample,nomatch = 0),2]
    if ( is.null(fontsize) ){
      fontsize = 10
    }
    if ( paired ){
      data2$patient <- condition[match(data2$sample,condition$sample,nomatch = 0),3]
      data2 <- data2[order(data2$Group,data2$patient),]
      p  <- ggboxplot(data2,x="Group",y="value",color = "Group",
                      id = "Type",palette = "jco",add = "jitter",facet.by= "Cell",
                      short.panel.labs = FALSE,line.color = "gray",line.size=0.4)+
        xlab("")+ylab("score")+
        geom_line(aes(group=patient),colour = "grey",linetype = 1,size = 0.2)+
        stat_compare_means(method = "wilcox.test",paired = T,label="p.format",label.x =1.5,label.y=1)+
        theme(axis.title.x = element_text(size=fontsize,face = "bold"))
    }else{
      p  <- ggboxplot(data2,x="Group",y="value",color = "Group",
                      id = "Type",palette = "jco",add = "jitter",facet.by= "Cell",
                      short.panel.labs = FALSE,line.color = "gray",line.size=0.4)+
        xlab("")+ylab("score")+
        stat_compare_means(method = "wilcox.test",label="p.format",label.x =1.5,label.y=1)+
        theme(axis.title.x = element_text(size=fontsize,face = "bold"))
    }
    return(p)
  }else{
    print(compare)
    stop("There are more than two groups")
  }
}
