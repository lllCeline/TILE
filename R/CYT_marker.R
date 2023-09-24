#' CYT marker
#' @description
#' The CYT (cytolytic activity) marker is estimated using the expression of granzyme A and perforin according to Rooney et al.,2015
#' @param exp_matrix a data matrix with rows labeling features gene symbol,and columns labeling samples,expression data should be normalized values,such as "TPM"
#' @examples
#' result <- CYT_calculation(demo_TPM)
#' @export

CYT_calculation <- function(exp_matrix){
  tpm <- data.frame(exp_matrix,check.names = FALSE)
  gene <- c("GZMA","PRF1")
  CYT_data <- tpm[which(rownames(tpm) %in% gene),]
  CYT_data <- as.data.frame(t(CYT_data))
  CYT_log2 <-log2(CYT_data+1)
  CYT_log2$CYT_value <- rowMeans(CYT_log2)
  CYT_log2 <- as.data.frame(t(CYT_log2))
  return(CYT_log2)
}
