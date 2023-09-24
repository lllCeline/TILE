#' IPRES marker
#' @description
#' Innately resistant tumors display a transcriptional signature (referred to as the IPRES or Innate anti-PD-1 Resistance) according to Hugo et al.,2016
#' @param exp_matrix a data matrix with rows labeling features gene symbol,and columns labeling samples
#' @param normalization a logical indicating whether expression data are normalized, such as "TPM",default is normalized values
#' @import GSVA
#' @examples
#' result <- IPRES(demo_TPM)
#' result <- IPRES(demo_counts,normalization = FALSE)
#' @export

IPRES <- function(exp_matrix,normalization = TRUE){
  input <- data.frame(exp_matrix,check.names = FALSE)
  gene_set <- ipres
  gene <- as.data.frame(t(gene_set[,2:ncol(gene_set)]),stringsAsFactors = F)
  colnames(gene) <- gene[1,];gene <- gene[-1,]
  ipres_gene <- c(gene)
  ipres_gene <- unlist(ipres_gene)
  ipres_gene <- as.data.frame(ipres_gene,stringsAsFactors = F)
  ipres_gene <- ipres_gene[which(ipres_gene[,1] != ""),]
  ipres_gene <- unique(ipres_gene)
  input$Gene <- rownames(input)
  input <- input[which(input$Gene %in% ipres_gene),]
  gene <- as.data.frame(table(input$Gene))
  gene <- gene[which(gene$Freq > 1),]
  if ( nrow(gene) >= 1 ){
    input_2 <- input[which( ! input$Gene %in% gene$Var1),]
    repeat_gene <- data.frame()
    for (i in 1:nrow(gene)){
      new <- input[which(input$Gene == gene[i,1]),]
      new <- as.data.frame(t(apply(new[,-1], 2, mean)))
      new$gene <- gene[i,1]
      if (i == 1){
        repeat_gene <- new
      }else{
        repeat_gene <- rbind(repeat_gene,new)
      }
    }
    input <- rbind(input_2,repeat_gene)
  }
  rownames(input) <- input$Gene
  input <- input[,-which(colnames(input)=="Gene")]
  input <- input[rowSums(input) > 0,]
  geneset <- split(gene_set[,3:ncol(gene_set)],gene_set[,2])
  rnull <- function(x){
    index <- which(x != "")
    x <- as.character(x[index])
  }
  geneset <- lapply(geneset, rnull)
  if (normalization){
    gsva_es <- gsva(log2(as.matrix(input)),geneset,method = "ssgsea",kcdf = "Gaussian")
  }else{
    gsva_es <- gsva(as.matrix(input),geneset,method = "ssgsea",kcdf = "Poisson")
  }
  return(gsva_es)
}
