#' novel_junction_ratios
#' @description Computes, per sample, the fraction of expressed splice junctions
#' that are novel (unannotated): x/(x+y), where x is the number of novel
#' junctions and y the number of annotated junctions, each counted only when
#' supported by more than `n` reads. With `global = TRUE` the counts are taken
#' across all junctions in the sample; with `global = FALSE` they are taken
#' within each gene.
#' @param JC junction count matrix, rows = "chr:start-end:strand", columns =
#' sample ids. A base matrix or a \pkg{Matrix} sparse matrix. Its rows must be
#' in the same order as `sjmeta` (same number of rows).
#' @param sjmeta output of `annotate_sj()`: the junction_meta data file
#' with added columns `gencode_gene_id`, `gencode_gene_name`.
#' Gene ids/names come from the packaged GENCODE release 29 models ([gencode_v29_genes]).
#' @param global (default=TRUE) set to TRUE for a global single score per sample,
#' FALSE for a computation on a per gene basis
#' @param n (default=5) a junction (novel or annotated) is counted only when it
#' has more than `n` reads
#' @return if `global` is TRUE, a data.frame with columns `sample_id` and `NJR`
#' (NA for samples with no annotated junction above `n` reads).
#' If FALSE, a numeric matrix with rows = genes ("gene_name:gene_id"),
#' columns = samples, entries = per-gene novel junction ratios. Only genes
#' with at least one annotated junction above `n` reads in some sample are
#' reported, and a gene's entry is NA in any sample where it has no annotated
#' junction above `n` reads.
#' @seealso [annotate_sj()], [gencode_v29_genes]
#' @importFrom Matrix sparseMatrix colSums
#' @export
novel_junction_ratios <- function(JC, sjmeta, global=TRUE, n=5){
  if(nrow(JC)!=nrow(sjmeta))
    stop('Error in novel_junction_ratios function. `JC` and `sjmeta` must have same number of rows')
  sj<-paste0(sjmeta$seqnames,':',sjmeta$start,'-',sjmeta$end,':',sjmeta$strand)
  if(!identical(sj,rownames(JC)))
    stop('Error in novel_junction_ratios function. `JC` and `sjmeta` must have corresponding splice junctions')

  # 0/1 indicator per junction x sample: does this junction count in this sample?
  qa<-which(sjmeta$annotated==1) #annotated junctions
  qu<-which(sjmeta$annotated==0) #novel junctions
  Ia<-JC[qa,,drop=FALSE]>n
  Iu<-JC[qu,,drop=FALSE]>n

  if(global){
    x<-colSums(Iu)
    y<-colSums(Ia)
    njr<-x/(x+y)
    njr[y==0]<-NA
    return(data.frame(sample_id=colnames(JC),NJR=as.numeric(njr)))
  }

  # gene x junction incidence matrix; junctions overlapping several genes
  # ("id1,id2") count towards each of them
  q<-which(!is.na(sjmeta$gencode_gene_id) & !is.na(sjmeta$gencode_gene_name))
  ids<-strsplit(sjmeta$gencode_gene_id[q],',')
  nms<-strsplit(sjmeta$gencode_gene_name[q],',')
  len<-lengths(ids)
  # junction metadata from older annotate_sj() versions de-duplicated the
  # name list (e.g. two RF00019 genes), so ids and names no longer pair up;
  # recover those names from the packaged GENCODE v29 gene models
  bad<-which(len!=lengths(nms))
  if(length(bad)>0){
    warning(length(bad),' junction(s) have gene id/name lists of different length; ',
            'taking their gene names from gencode_v29_genes')
    gv29<-prepTCGAdata::gencode_v29_genes
    nms[bad]<-lapply(ids[bad],function(v) gv29$gene_name[match(v,gv29$gene_id)])
  }
  genes<-paste(unlist(nms),unlist(ids),sep=':')
  junc<-rep(q,len)
  ug<-unique(genes)
  G<-sparseMatrix(i=match(genes,ug),j=junc,x=1,
                  dims=c(length(ug),nrow(JC)),dimnames=list(ug,NULL))
  G<-(G>0)*1  # guard against a gene listed twice for the same junction

  y<-as.matrix(G[,qa,drop=FALSE] %*% Ia)
  # report only genes with an annotated junction above n reads in some sample
  keep<-which(rowSums(y)>0)
  y<-y[keep,,drop=FALSE]
  x<-as.matrix(G[keep,qu,drop=FALSE] %*% Iu)
  out<-x/(x+y)
  out[y==0]<-NA
  colnames(out)<-colnames(JC)
  return(out)
}
