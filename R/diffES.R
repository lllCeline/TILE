#' limma analysis of cell type score
#' @description
#' Different immune cell types were analyzed by the limma package
#' @param score_matrix a data matrix with rows labeling features immune cell type, and columns labeling samples.
#' @param groupinfo a matrix with two columns: 'sample' column labeling sample id,'groups' column labeling group information
#' @import limma
#' @examples
#' result <- diffES(score_matrix = cell_enrich,groupinfo = demo_groupinfo)
#' @export


diffES <- function(score_matrix,groupinfo){
  gsva_es <- data.frame(score_matrix,check.names = FALSE)
  condition <- data.frame(groupinfo,check.names = FALSE)
  compare <- unique(condition$groups)
  if (TRUE %in% (condition$sample %in% rownames(gsva_es))){
    gsva_es <- data.frame(t(gsva_es),check.names = FALSE)
  }
  if (length(compare) == 2){
    group1 <- compare[1]
    group2 <- compare[2]
    gsva_es <- gsva_es[,which(colnames(gsva_es) %in% condition$sample)]
    f<- factor(condition$groups,levels = compare)
    design <- model.matrix(~0+f)
    colnames(design) <- levels(f)
    rownames(design) <- condition$sample
    fit <- lmFit(gsva_es,design)
    if (all(fit$df.residual) == 0){
      res <- as.data.frame(gsva_es)
      res$log2FoldChange <- log2(res[,1]/res[,2])
      return(res)
    }else{
      contmatrix <- makeContrasts(contrasts = paste(levels(f)[1],levels(f)[2],sep = "-"),levels = design)
      fit <- contrasts.fit(fit,contmatrix)
      fit <- eBayes(fit = fit)
      result <- topTable(fit,coef = 1,adjust.method  = "BH",sort.by = "p" ,number = nrow(gsva_es))
      output <- data.frame(gsva_es[rownames(result),],result,stringsAsFactors = F,check.names = F,row.names = rownames(gsva_es))
      output <- output[,-c(ncol(output),ncol(output)-3,ncol(output)-4)]
      colnames(output)[(ncol(output)-2):ncol(output)] <- c("log2FoldChange","pvalue","padj")
      return(output)
    }
  }else{
    print(compare)
    stop("There are more than two groups")
  }
}
