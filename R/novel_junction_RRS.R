#' novel_junction_RRS
#' @description Computes the relative read support of novel junctions within a gene,
#'  X/(X+Y), where X=number of reads for the novel junction, 
#'  Y is the median read support for annotated junctions within the gene.
#' @param JC junction count matrix, rows = "chr:start-end:strand", columns =
#' sample ids. A base matrix or a \pkg{Matrix} sparse matrix. Its rows must be
#' in the same order as `sjmeta` (same number of rows).
#' @param sjmeta output of `annotate_sj()`: the junction_meta data file 
#' with added columns `gencode_gene_id`, `gencode_gene_name`. 
#' Gene ids/names come from the packaged GENCODE release 29 models ([gencode_v29_genes]).
#' @return a matrix in same format as JC, but subset to just the novel junctions,
#' with entries showing the RRS scores of these junctions
#' @seealso [annotate_sj()], [gencode_v29_genes]
#' @importFrom stats median
#' @importFrom tidyr separate_rows
#' @export
novel_junction_RRS <- function(JC, sjmeta){
  
  #--- compute RRS scores --------------------------------------------------#
  sjmeta$sj<-rownames(JC)
  
  q<-which(!is.na(sjmeta$gencode_gene_name))
  sjmeta<-sjmeta[q,]
  JC<-JC[q,]

  qa<-which(sjmeta$annotated==1) #annotated junctions
  JC.a<-JC[qa,]
  sjmeta.a<-sjmeta[qa,]

  q<-which(sjmeta$annotated==0) #novel junctions
  JC.n<-JC[q,]
  sjmeta.n<-sjmeta[q,]
  out<-0*JC.n
  
  start_time <- Sys.time()
  # A carriage-return progress bar (txtProgressBar) only redraws in a live
  # terminal; it shows nothing useful once output is redirected to a log file,
  # run as an RStudio background Job, or captured by another process (all
  # common for a run this long). Report progress as plain, whole lines instead
  # - message() always writes a complete, immediately-flushed line, so it's
  # visible wherever the output ends up. ~20 updates regardless of ncol(JC).
  report_every <- max(1L, round(ncol(JC.n) / 20))
  for(i in 1:ncol(JC.n)){
    k<-which(JC.n[,i]>0)
    df.x<-dplyr::select(sjmeta.n[k,],c(gencode_gene_id,sj))
    df.x$count<-JC.n[k,i]
    df.x <- separate_rows(df.x, gencode_gene_id, sep = ",")
    
    k<-which(JC.a[,i]>0)
    df.y<-dplyr::select(sjmeta.a[k,],c(gencode_gene_id,sj))
    df.y$count<-JC.a[k,i]
    df.y <- tidyr::separate_rows(df.y, gencode_gene_id, sep = ",")
    df.y <- filter(df.y,count>1)
    
    m<-df.y %>% 
      group_by(gencode_gene_id) %>% 
      summarize(med=max(c(1,median(count))))
    
    z<-intersect(df.x$gencode_gene_id,df.y$gencode_gene_id)
    xy<-left_join(filter(df.x,is.element(gencode_gene_id,z)),
                  filter(m,is.element(gencode_gene_id,z)),
                  by='gencode_gene_id')
    xy$RRS<-xy$count/(xy$count+xy$med)
    out[xy$sj,i]<-xy$RRS
    if (i %% report_every == 0 || i == ncol(JC.n))
      message(sprintf("  novel_junction_RRS: sample %d/%d (%.0f%%)", i, ncol(JC.n), 100 * i / ncol(JC.n)))
  }
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  message(sprintf("novel_junction_RRS: completed %d sample(s) in %.2f seconds.", ncol(JC.n), elapsed))
  return(out)#list(RRS=out,RRSmeta=sjmeta.n))
}



