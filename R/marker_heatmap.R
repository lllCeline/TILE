#' heatmap visualization
#' @description
#' supervised analyses of cell type score between groups
#' @param score_matrix a data matrix with rows labeling features immune cell type,and columns labeling samples
#' @param groupinfo a matrix with two columns:'sample' column labeling sampleid,'groups' column labeling group information
#' @param scaled a logical indicating whether the values would be centered and scaled in either the row direction or the column direction, or none
#' @param title_size a numeric value, change the size of title, default is NULL(function can provide the font size of title automatically)
#' @param rowsize a numeric value, change the size of row labels, default is NULL(function can provide the font size of row labels automatically)
#' @param colsize a numeric value, change the size of column labels, default is NULL(function can provide the font size of column labels automatically)
#' @import ComplexHeatmap
#' @import ggpubr
#' @import grid
#' @examples
#' p <- marker_heatmap(score_matrix = cell_enrich,groupinfo = demo_groupinfo, title_size = 8)
#' print(p)
#' @export

marker_heatmap <- function(score_matrix,groupinfo,scaled = "row",title_size = NULL,colsize = NULL,rowsize = NULL){
  input <- data.frame(score_matrix, check.names = FALSE)
  condition <- data.frame(groupinfo,check.names = FALSE)
  if (TRUE %in% (condition$sample %in% colnames(input))){
    input <- data.frame(t(input),check.names = FALSE)
  }
  scale_method <- c("row","col","none")
  if (!scaled %in% scale_method){
    stop('"scaled" should be one of "row","col","none"')
  }
  if (scaled == "row"){
    tmp <- apply(input, 2, scale)
    tmp <- data.frame(tmp,row.names = rownames(input),check.names = FALSE)
    input <- tmp
  }else if (scaled == "col"){
    tmp <- apply(input,1,scale)
    tmp <- data.frame(tmp,row.names = rownames(input),check.names = FALSE)
    input <- tmp
  }
  cell_type <- colnames(input)
  input <- input[which(rownames(input) %in% condition$sample),]
  input$group <- condition[match(rownames(input),condition$sample,nomatch = 0),2]
  input <- input[order(input$group),]
  group <- unique(condition$groups)
  color_n <- length(group)
  if(color_n==1){
    col_nw = "darkorange"
  }else if (color_n == 2){
    col_nw =c("#CC0033","#003399")
  }else{col_nw = RColorBrewer::brewer.pal(n=color_n, name="Set3")}
  for ( i in 1:length(group)){
    data_sub <- input[which(input$group == group[i]),]
    data_sub <- as.matrix(data_sub[,-ncol(data_sub),drop = FALSE])
    h <- Heatmap(data_sub,cluster_columns = FALSE);h <- draw(h)
    order <- as.matrix(row_order(h))
    assign(paste(group[i],"_sample",sep = ""),rownames(data_sub)[order])
  }
  input <- as.matrix(t(input[,-ncol(input)]))
  horder <- c()
  for ( i in 1:length(group)){
    sample <- get(paste0(group[i],"_sample"))
    horder <- append(horder,sample)
  }
  group.assign <- setNames(col_nw,group)
  match.id <- match(colnames(input),condition$sample,nomatch = 0)
  df <- condition[match.id,2,drop = FALSE]
  if (is.null(colsize)){
    if ( ncol(input) <= 15 ){
      colsize = 10
    }else if (ncol(input) > 15 & ncol(input) <= 30){
      colsize = 8
    }else{colsize = 6}
  }
  if (is.null(title_size)){
    if ( ncol(input) <= 15 ){
      title_size = 20
    }else if (ncol(input) > 15 & ncol(input) <= 30){
      title_size = 22
    }else{title_size = 28}
  }
  if (is.null(rowsize)){
    if ( nrow(input) <= 15 ){
      rowsize = 10
    }else if (nrow(input) > 15 & nrow(input) <= 30){
      rowsize = 8
    }else{rowsize = 6}
  }
  legend_size = max(rowsize,colsize)-2
  ha <- HeatmapAnnotation(df=df,
                          col = list(groups=group.assign),
                          annotation_height = unit(0.5,"cm"),
                          annotation_label = "Group",
                          show_annotation_name = FALSE,
                          annotation_legend_param = list(title_gp = gpar(fontsize = legend_size),labels_gp = gpar(fontsize = legend_size)))
  p <- Heatmap(input,
               name = " ",
               column_title = "cell type heatmap",
               column_title_gp = gpar(fontsize = title_size,fontface = "bold",col = "black"),
               heatmap_legend_param = list(title = "",legend_height = unit(3, "cm"),labels_gp = gpar(fontsize = legend_size)),
               row_names_gp = gpar(fontsize = rowsize),
               column_names_gp = gpar(fontsize = colsize),
               cluster_columns = FALSE,
               bottom_annotation = ha,
               column_order = horder ,
               row_names_max_width = max_text_width(rownames(input),gp = gpar(fontsize = rowsize)))
  return(p)
}
