#' box plot visualization
#' @description
#' Box plots for immnume signatures score with clinical features
#' @param score_matrix a data matrix with rows labeling features immune cell type,and columns labeling samples
#' @param paired a logical indicating whether you want a paired wilcoxon test, which is used to compare the means of two samples when each observation in one sample can be paired with an observation in the other sample,default is "FALSE"
#' @param groupinfo a matrix with two columns:'sample' column labeling sampleid,'groups' column labeling group information.
#' when @param paired was setted as'TRUE',@param groupinfo need to provide three columns,the added column 'patient' labeling the paired relation of samples
#' @param fontsize a numeric value, change the font size, default is NULL(function can provide the font size automatically)
#' @import ggplot2
#' @import ggpubr
#' @import reshape2
#' @examples
#' p <- marker_boxplot(score_matrix = cell_enrich,groupinfo = demo_groupinfo,fontsize = 10,paired = FALSE)
#' print(p)
#' #paired sample
#' p <- marker_boxplot(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
#' print(p)
#' @export

marker_boxplot <- function(score_matrix,groupinfo,paired = FALSE,fontsize = NULL){
  data <- data.frame(score_matrix,check.names = FALSE)
  condition <- data.frame(groupinfo,check.names = FALSE)
  if (TRUE %in% (condition$sample %in% colnames(data))){
    data <- t(data)
  }
  cell_type <- colnames(data)
  data <- data[which(rownames(data) %in% condition$sample),]
  data <- as.data.frame(t(data))
  compare <- unique(condition$groups)
  if (length(compare) == 2){
    group1 <- compare[1]
    group2 <- compare[2]
    condition <- condition[which(condition$groups == group1|condition$groups == group2),]
    data <- subset(data,select = as.character(condition$sample))
    data$gene <- rownames(data)
    data2 <- reshape2::melt(data,variable.name = "sample")
    data2$Group <- condition[match(data2$sample,condition$sample,nomatch = 0),2]
    if (is.null(fontsize)){
      if (ncol(data2) <= 30){
        fontsize = 10
      }else{
        fontsize = 8
      }
    }
    if (paired){
      data2$patient <- condition[match(data2$sample,condition$sample,nomatch = 0),3]
      data2 <- data2[order(data2$Group,data2$patient),]
      p <- ggplot(data2,aes(x=gene,y=value,fill=Group))+
        geom_boxplot()+
        stat_boxplot(geom = "errorbar")+
        theme_bw()+
        theme(axis.text.x = element_text(size = fontsize,angle = 90,vjust = 0.5,hjust = 1))+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        labs(y="score",x="")+
        theme(plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values = c("#EF0808","#085A9C"))+
        stat_compare_means(method = "wilcox.test",paired = TRUE,label="p.signif")
    }else{
      p <- ggplot(data2,aes(x=gene,y=value,fill=Group))+
        geom_boxplot()+
        stat_boxplot(geom = "errorbar")+
        theme_bw()+
        theme(axis.text.x = element_text(size = fontsize,angle = 90,vjust = 0.5,hjust = 1))+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        labs(y="score",x="")+
        theme(plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values = c("#EF0808","#085A9C"))+
        stat_compare_means(method = "wilcox.test",label="p.signif")
    }
  }else{
    print(compare)
    stop("There are more than two groups")
  }
  return(p)
}
