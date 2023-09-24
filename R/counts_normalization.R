#' raw reads count normalization
#' @param counts_matrix a data matrix with rows labeling features gene symbol,and columns labeling samples,expression data shoule be raw reads count values
#' @param method a character string provide the normalization methods, can be one of c("TPM","FPKM")
#' @param gene_len a character string provide the gene length data,gene type is SYMBOL,gene length database can selected from c("hg19_gene_symbol_length","hg38_gene_symbol_length")
#' @examples
#' result <- counts_normlization(demo_counts,method = "FPKM",gene_len = "hg38_gene_symbol_length")
#' @export

counts_normlization <- function(counts_matrix,method = "TPM",gene_len = "hg19_gene_symbol_length"){
  input <- data.frame(counts_matrix,check.names = FALSE)
  method_list <- c("TPM","FPKM")
  if (!method %in% method_list){
    stop('"method" should be one of "TPM","FPKM".')
  }
  gene_len_list <- c("hg19_gene_symbol_length","hg38_gene_symbol_length")
  if (!gene_len %in% gene_len_list){
    stop('"gene_len" should be one of "hg19_gene_symbol_length","hg38_gene_symbol_length".')
  }
  glen <- get(gene_len)
  glen <- data.frame(glen,check.names = FALSE)
  input <- input[which(rownames(input) %in% glen$gene_symbol),]
  if (method == "TPM"){
    tpm <- as.data.frame(matrix(nrow=nrow(input),ncol=0),row.names = rownames(input))
    for (i in 1:ncol(input)){
      x <- input[,i,drop = FALSE]
      y <- data.frame(x,glen[match(rownames(input),glen$gene_symbol,nomatch = 0),2])
      y$rate <- apply(y,1,function (j) j[1]/j[2])
      total <- apply(y,2,sum)
      tpm[,i] <-  data.frame(round(y$rate/total["rate"]*1e6,3))
      colnames(tpm)[i] <- colnames(input)[i]
    }
    return(tpm)
  }
  else if (method == "FPKM") {
    fpkm <- as.data.frame(matrix(nrow=nrow(input),ncol=0),row.names = rownames(input))
    for (i in 1:ncol(input)){
      x <- input[,i,drop = FALSE]
      total <- lapply(x, sum)
      total <- as.numeric(total)
      rate <- 1e9 / total
      y <- data.frame(x,glen[match(rownames(input),glen$gene_symbol,nomatch = 0),2])
      fpkm[,i] <- as.data.frame(round(apply(y, 1, function(j) j[1] * rate/j[2]),3))
      colnames(fpkm)[i] <- colnames(input)[i]
    }
    return(fpkm)
  }
}
