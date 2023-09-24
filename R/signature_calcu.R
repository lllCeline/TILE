#' signature calculation
#' @description
#' This function can be used to calculate immnune signature
#' @param exp_matrix a data matrix with rows labeling features gene symbol,and columns labeling samples,expression data should be normalized values,such as "TPM"
#' @param db a character string can be 'integration_signature' or 'tls_sig'
#' integration_signature integrate tumor microenvironment-related signature,evaluated the abundance of immune cell populations in the tumor microenvironment according to Pornpimol et al.,2017 and Ott et al.,2019
#'
#' tls_sig signature can be used to calculate the score of Tertiary lymphoid structures according to Tokunaga et al.,2020 and Rita et al.,2020
#' @examples
#' result <- sig_calculate(exp_matrix = demo_TPM,db = "tls_sig")
#' @export


sig_calculate <- function(exp_matrix,db){
  db_list = c("integration_signature","tls_sig")
  if (!db %in% db_list){
    stop('"cell_type" should be one of "integration_signature","tls_sig".')
  }
  tpm <- data.frame(exp_matrix,check.names = FALSE)
  geneset <- get(db)
  data <- tpm[which(rownames(tpm) %in% geneset$Gene),]
  data_log2 <- as.data.frame(log2(data+1))
  data_log2$Gene <- rownames(data_log2)
  data_log2$cell_type <- geneset[match(data_log2$Gene,geneset$Gene),2]
  data_log2 <- distinct(data_log2)
  score <- data.frame()
  for ( i in 1:(ncol(data_log2)-2)){
    sample_name <- colnames(data_log2)[i]
    cell_type_score <- group_by(data_log2[,c(i,ncol(data_log2))],cell_type)
    colnames(cell_type_score)[1] <- "sample"
    cell_type_score <-  cell_type_score %>%
      summarise(sample = mean(sample)) %>%
      data.frame()
    colnames(cell_type_score)[2] <- sample_name
    if (i == 1){
      score <- cell_type_score
    }else{
      score <- merge(score,cell_type_score,by = "cell_type")
    }
  }
  rownames(score) <- cell_type_score$cell_type
  score <- score[,-1]
  return(score)
}
