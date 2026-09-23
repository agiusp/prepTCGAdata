#' count_novel_sjs_per_pathway
#' @description Rolls [count_novel_sjs()]'s per-gene novel-junction counts up
#' to per-pathway (gene-set) scores, per sample: for each pathway, the sum of
#' its member genes' novel-junction counts.
#' @param nsj output of [count_novel_sjs()]: a gene x sample matrix whose row
#' names are `"gene_name:gene_id"`.
#' @param pways gene sets to score, as a named list - each element a character
#' vector of gene names or gene ids (matched against `nsj`'s `gene_name`), the
#' list names giving the pathway names.
#' @return a matrix, one row per pathway in `pways` and one column per
#' `sample_id` in `nsj`, with entries showing the summed novel-junction count
#' of that pathway's member genes, for that sample. A pathway with no member
#' gene present in `nsj` still gets a row (all zero).
#' @seealso [get_tcga_data()]
#' @importFrom Matrix sparseMatrix t
#' @export
count_novel_sjs_per_pathway <- function(nsj, pways){
  
  pathway_sums <- function(expr, df, sig_list) {
    
    genes <- df$gene_name
    n_genes <- length(genes)
    n_pathways <- length(sig_list)
    
    # Build sparse gene x pathway indicator matrix
    # (avoids allocating a dense 50000 x n_pathways matrix)
    gene_idx <- match(unlist(sig_list), genes)
    pathway_idx <- rep(seq_along(sig_list), lengths(sig_list))
    
    # Drop genes in signatures that aren't in your expression matrix
    keep <- !is.na(gene_idx)
    gene_idx <- gene_idx[keep]
    pathway_idx <- pathway_idx[keep]
    
    indicator <- sparseMatrix(
      i = gene_idx,
      j = pathway_idx,
      x = 1,
      dims = c(n_genes, n_pathways)
    )
    colnames(indicator) <- names(sig_list)
    
    # One matrix multiplication does all the summing:
    # (pathways x genes) %*% (genes x samples) = pathways x samples
    pathway_scores <- as.matrix(t(indicator) %*% expr)
    
    pathway_scores
  }
  
  g<-unlist(strsplit(rownames(nsj),':'))
  df<-data.frame(gene_name=g[seq(1,length(g),2)],
                 gene_id=g[seq(2,length(g),2)])
  return(pathway_sums(nsj,df,pways))
}
