#' TILs marker
#' @description
#' Total TILs score and cell type enrichment scores are calculated according to Danaher et al.,2017
#' @param exp_matrix a data matrix with rows labeling features gene symbol, and columns labeling samples,expression data should be normalized values,such as "TPM"
#' @import dplyr
#' @examples
#' result <- Co_expression(demo_TPM)
#' TILs_score <- data.frame(result[[1]],check.names = FALSE)
#' cell_enrich <- data.frame(result[[2]],check.names = FALSE)
#' @export

Co_expression <- function(exp_matrix){
  exp <- data.frame(exp_matrix,check.names = FALSE)
  cell_type <- co_expression_selected_gene
  marker <- exp[which(rownames(exp) %in% cell_type$Gene),]
  marker_nrom <- log2(marker)
  marker_nrom <- apply(marker_nrom, 1, function(x){x <- ifelse (x == "-Inf",0,x)})
  marker_nrom <- as.data.frame(t(marker_nrom))
  marker_nrom$cell_type = cell_type[match(rownames(marker_nrom),cell_type$Gene,nomatch = 0),2]
  cell_type_score_matrix <- data.frame(check.names = F)
  for (i in 1:(ncol(marker_nrom)-1)){
    name <- colnames(marker_nrom)[i]
    cell_type_score <- group_by(marker_nrom[,c(i,ncol(marker_nrom))],cell_type)
    colnames(cell_type_score) <- c("sample","cell_type")
    cell_type_score <-  cell_type_score %>%
      summarise(sample = mean(sample)) %>%
      data.frame()
    colnames(cell_type_score)[2] <- name
    if (i == 1){
      cell_type_score_matrix <- cell_type_score
    }else{
      cell_type_score_matrix <- merge(cell_type_score_matrix,cell_type_score,by = "cell_type")
    }
  }
  rownames(cell_type_score_matrix) <- cell_type_score_matrix$cell_type
  TIL_matrix <- cell_type_score_matrix[,-1]
  TIL_matrix <- TIL_matrix[cor(t(TIL_matrix))[,"CD45"]>0.6,]
  TIL_matrix["TILs_Score",] <- apply(TIL_matrix, 2, mean)
  TIL_matrix <- rbind(cell_type_score_matrix[,-1],TIL_matrix["TILs_Score",])
  diff_data <- TIL_matrix
  TIL_matrix <- data.frame(t(TIL_matrix),check.names = FALSE)
  cell_enrich <- TIL_matrix*0
  for (i in 1:ncol(cell_enrich)){
    if (colnames(cell_enrich)[i] == "TILs_Score"){
      cell_enrich[,i] <- TIL_matrix[,i]
    }else{
      mod <- lm(TIL_matrix[,i]~TIL_matrix[,"TILs_Score"])$coef
      cell_enrich[,i] <- TIL_matrix[,i]-(mod[1]+mod[2]*TIL_matrix[,"TILs_Score"])
    }
  }
  cell_enrich <- t(cell_enrich)
  outlist <- list(TIL_score = diff_data,cell_enrich = cell_enrich)
  return(outlist)
}
