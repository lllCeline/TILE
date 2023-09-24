#' GEP marker
#' @description
#' T-cell–inflamed gene-expression profile (GEP),which is a biomarker that can predict response to anti–programmed cell death 1 (PD-1) therapy across multiple tumor types according to Ayers, Mark et al.,2017
#' @param exp_matrix a data matrix with rows labeling features gene symbol,and columns labeling samples
#' @import dplyr
#' @examples
#' result <- GEP(exp_matrix = demo_TPM)
#' @export

GEP <- function(exp_matrix){
  input <- data.frame(exp_matrix,check.names = FALSE)
  input <- log10(input+1)
  if (TRUE %in% (gep_gene$gene %in% colnames(input))){
    input <- data.frame(t(input),check.names = FALSE)
  }
  gep_tpm <- input[which(rownames(input) %in% gep_gene$gene),]
  hk <- gep_hk
  hk_tpm <- input[which(rownames(input) %in% hk$gene),]
  all <- as.data.frame(t(rbind(gep_tpm,hk_tpm)))
  all_hk_aver <- all %>%
    mutate(aver = rowMeans(dplyr::select(.,hk$gene)))
  all_hk_nrom <- sweep(all_hk_aver,1,all_hk_aver$aver,"-")
  gep_tpm_nrom <- all_hk_nrom[,which(colnames(all_hk_nrom) %in% gep_gene$gene)]
  gep_gene$weight <- as.numeric(gep_gene$weight)
  calu <- function(x){
    s <- 0
    for (i in 1:ncol(gep_tpm_nrom)){
      i = 1
      c <- gep_gene[which(colnames(gep_tpm_nrom) %in% colnames(gep_tpm_nrom)[1]),2] * x[i]
      s <- s + c
    }
    return(s)
  }
  score <- apply(gep_tpm_nrom,1,calu)
  gep_score <- data.frame(GEP_score = score)
  gep_score <- as.data.frame(t(gep_score))
  return(gep_score)
}
