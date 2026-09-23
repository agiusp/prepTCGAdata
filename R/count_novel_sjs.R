#' count_novel_sjs
#' @description Counts, per gene and per sample, the novel (unannotated) splice
#' junctions with non-zero reads.
#' @param JC junction count matrix, rows = "chr:start-end:strand", columns =
#' sample ids. A base matrix or a \pkg{Matrix} sparse matrix. Its rows must be
#' in the same order as `sjmeta` (same number of rows).
#' @param sjmeta output of `annotate_sj()`: a data.frame with columns
#' `gencode_gene_id`, `gencode_gene_name` and an annotation-status column
#' (`annotated`). Gene ids/names come from the packaged GENCODE release 29
#' models ([gencode_v29_genes]).
#' @return a \pkg{Matrix} object, one row per gene (row names `"gene_name:gene_id"`)
#' and one column per `sample_id` in `JC`, with entries showing the count of
#' novel (unannotated) splice junctions with non-zero reads intersecting that
#' gene, for that sample. Genes with a zero count in every sample are dropped.
#' @seealso [annotate_sj()], [gencode_v29_genes]
#' @importFrom Matrix Matrix
#' @export
count_novel_sjs <- function(JC, sjmeta){
  q<-intersect(which(!is.na(sjmeta$gencode_gene_id)),
               which(!is.na(sjmeta$gencode_gene_name)))
  gannot<-data.frame(gene_name=unlist(strsplit(sjmeta$gencode_gene_name[q],',')),
                     gene_id=unlist(strsplit(sjmeta$gencode_gene_id[q],',')))
  gannot<-unique(gannot)
  rownames(gannot)<-gannot$gene_id
  out<-Matrix(0,nrow(gannot),ncol(JC))
  rownames(out)<-paste(gannot$gene_name,gannot$gene_id,sep=':')
  colnames(out)<-colnames(JC)
  
  start_time <- Sys.time()
  # A carriage-return progress bar (txtProgressBar) only redraws in a live
  # terminal; it shows nothing useful once output is redirected to a log file,
  # run as an RStudio background Job, or captured by another process (all
  # common for a run this long). Report progress as plain, whole lines instead
  # - message() always writes a complete, immediately-flushed line, so it's
  # visible wherever the output ends up. ~20 updates regardless of ncol(JC).
  report_every <- max(1L, round(ncol(JC) / 20))

  for(i in 1:ncol(JC)){
    df<-filter(sjmeta[intersect(q,which(JC[,i]>0)),],annotated==0)
    tmp<-data.frame(id=unlist(strsplit(df$gencode_gene_id,',')),
                    name=unlist(strsplit(df$gencode_gene_name,','))) %>%
      group_by(id,name) %>% summarize(count=n(),.groups = "drop_last")
    tmp$rname=paste(tmp$name,tmp$id,sep=':')
    out[tmp$rname,i]<-tmp$count
    if (i %% report_every == 0 || i == ncol(JC))
      message(sprintf("  count_novel_sjs: sample %d/%d (%.0f%%)", i, ncol(JC), 100 * i / ncol(JC)))
  }
  n<-apply(out,1,sum)
  out<-out[which(n>0),]

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  message(sprintf("count_novel_sjs: completed %d sample(s) in %.2f seconds.", ncol(JC), elapsed))
  return(out)
}
  
#   split_vector_to_matrix <- function(vec) {
#     # Split each string by comma
#     split_list <- strsplit(vec, ",")
#     
#     # Find the max number of elements in any split string
#     max_cols <- max(sapply(split_list, length))
#     
#     # Pad each vector with NA up to max_cols, then rbind
#     padded <- lapply(split_list, function(x) {
#       length(x) <- max_cols  # extends with NA automatically
#       x
#     })
#     
#     # Combine into a matrix (rows = original elements)
#     result <- do.call(rbind, padded)
#     rownames(result) <- NULL
#     
#     return(result)
#   }
#   
#   #--- build gene USJ counts --------------------------------------------------#
#   
#   q<-which(!is.na(sjmeta$gencode_gene_name))
#   q<-intersect(q,which(sjmeta$annotated==0)) #novel junctions
#   JC<-JC[q,]
#   sjmeta<-sjmeta[q,]
#   q.comma<-grep(',',sjmeta$gencode_gene_id)
#   out<-NULL
#   for(i in 1:ncol(JC)){
#     k<-which(JC[,i]>0)
#     k1<-setdiff(k,q.comma)
#     tt1<-sjmeta[k1,] %>% 
#       group_by(gencode_gene_id,gencode_gene_name) %>% 
#       summarize(count=n(), .groups = "drop")
#     k2<-setdiff(k,k1)
#     ggi<-split_vector_to_matrix(sjmeta$gencode_gene_id[k2])
#     ggn<-split_vector_to_matrix(sjmeta$gencode_gene_name[k2])
#     tt<-tt1
#     for(j in 1:ncol(ggn)){
#       tt.tmp<-data.frame(gencode_gene_id=ggi[,j],
#                          gencode_gene_name=ggn[,j])
#       tt2<-tt.tmp %>% 
#         group_by(gencode_gene_id,gencode_gene_name) %>% 
#         summarize(count2=n(), .groups = "drop")
#       tt<-full_join(tt,tt2,by=c('gencode_gene_id','gencode_gene_name'))
#       tt$count[is.na(tt$count)] <- 0
#       tt$count2[is.na(tt$count2)] <- 0
#       tt$count<-tt$count+tt$count2
#       tt<-select(tt,!count2)
#     }
#     colnames(tt)[3]<-colnames(JC)[i]
#     if(i==1){
#       out<-tt
#     } else {
#       out<-full_join(out,tt,by=c('gencode_gene_id','gencode_gene_name'))
#     }
#     pb <- txtProgressBar(min = 1, max = ncol(JC), style = 3)
#   }
#   
#   return(out)
# }


