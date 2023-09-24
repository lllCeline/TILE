#' ssGESA method
#' @description
#' use single sample gene set enrichment analysis (ssGSEA) (Barbie et al., 2009) to identify immune cell types that are over-represented in the tumor microenvironment
#' @param exp_matrix a data matrix with rows labeling features gene symbol,and columns labeling samples
#' @param normalization a logical indicating whether expression data are normalized, such as "TPM",default is normalized values
#' @param db a character string can be one of c("msigdb_C2","msigdb_C5","msigdb_C6","msigdb_C7","msigdb_H","cancerSEA","TILs_28_subpopulations")
#' msigdb db:https://www.gsea-msigdb.org/gsea/msigdb
#' cancerSEA:http://biocc.hrbmu.edu.cn/CancerSEA
#' TILs_28_subpopulations:pan-cancer immune cell subpopulations according to Pornpimol et al.,2017
#' @import GSVA
#' @examples
#' result <- ssgsea(demo_counts,normalization = FALSE,db = "cancerSEA")
#' result <- ssgsea(demo_TPM,normalization = TRUE,db = "msigdb_H")
#' @export

ssgsea <- function(exp_matrix,normalization = TRUE,db){
  db_list = c("msigdb_C2","msigdb_C5","msigdb_C6","msigdb_C7","msigdb_H","cancerSEA","TILs_28_subpopulations")
  if (!db %in% db_list){
    stop('"db" should be one of "msigdb_C2","msigdb_C5","msigdb_C6","msigdb_C7","msigdb_H","cancerSEA","TILs_28_subpopulations".')
  }
  #input <- data.frame(exp_matrix,check.names = FALSE)
  input <- as.matrix(exp_matrix)
  if (grepl("msigdb",db)){
    geneset <- get(db)
  }else{
    geneset <- get(db)
    geneset <- split(geneset$Metagene,geneset$Cell_type)
  }
  if (normalization){
    input <- log2(input+1)
    gsva_es <- gsva(input,geneset,method="ssgsea",kcdf = "Gaussian")
  }else{
    gsva_es <- gsva(input,geneset,method="ssgsea",kcdf = "Poisson")
  }
  return(gsva_es)
}
