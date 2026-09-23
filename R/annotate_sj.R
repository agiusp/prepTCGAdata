#' annotate_sj
#' @description Annotates a table of splice-junction (or any) genomic intervals
#' with the GENCODE gene ids and gene names whose gene body overlaps each
#' interval. Two columns, `gencode_gene_id` and `gencode_gene_name`, are added
#' to `sjmeta`; when an interval overlaps several genes the ids/names are
#' collapsed into a single comma-separated string (NA when nothing overlaps).
#' The two strings are kept strictly parallel: `gencode_gene_id` and
#' `gencode_gene_name` always have the same number of comma-separated entries
#' and the i-th name is the name of the i-th id. A gene symbol shared by several
#' distinct loci (e.g. the Rfam small RNAs such as `RF00019`) therefore appears
#' once per overlapping gene id.
#'
#' The GENCODE reference is the release 29 gene set shipped with the package
#' ([gencode_v29_genes]); there is no GENCODE file to download or point at.
#' Release 29 is used because recount3 builds its TCGA splice junctions against
#' `"gencode_v29"`.
#' @param sjmeta the junction metadata: either a path/URL (a delimited text
#' file, or a serialized data.frame `.rds`), or an already-loaded data.frame /
#' list. Either way it needs (at least) the columns `start`, `end` and a
#' chromosome column (`chr`, or `seqnames` if `chr` is absent). `start`/`end`
#' are treated as 1-based inclusive coordinates and all original columns are
#' retained.
#' @return a data.frame: `sjmeta` with `gencode_gene_id` and
#' `gencode_gene_name` appended.
#' @seealso [gencode_v29_genes]
#' @export
annotate_sj <- function(sjmeta) {

  # ---- 1. read the interval table -----------------------------------------
  if (is.data.frame(sjmeta) || is.list(sjmeta)) {
    # already loaded (e.g. the caller did `readRDS()` themselves)
    sj <- as.data.frame(sjmeta)
  } else {
    if (!is.character(sjmeta) || length(sjmeta) != 1L || is.na(sjmeta)) {
      stop(
        "sjmeta must be a single file path/URL (character(1)), or an ",
        "already-loaded data.frame - got ", class(sjmeta)[1],
        if (is.character(sjmeta)) paste0(" of length ", length(sjmeta)) else ""
      )
    }
    # .rds -> serialized data.frame; otherwise a delimited text file / URL.
    is_rds <- grepl("\\.rds$", sjmeta, ignore.case = TRUE)
    if (is_rds) {
      con <- if (grepl("^(https?|ftp)://", sjmeta)) gzcon(url(sjmeta)) else sjmeta
      sj  <- as.data.frame(readRDS(con))
      if (inherits(con, "connection")) close(con)
    } else {
      sj <- as.data.frame(data.table::fread(sjmeta))
    }
  }

  chr_col <- intersect(c("chr", "seqnames"), colnames(sj))[1]
  missing <- setdiff(c("start", "end"), colnames(sj))
  if (is.na(chr_col)) missing <- c("chr/seqnames", missing)
  if (length(missing))
    stop("sjmeta is missing required column(s): ", paste(missing, collapse = ", "))

  sj_gr <- GenomicRanges::GRanges(
    seqnames = as.character(sj[[chr_col]]),
    ranges   = IRanges::IRanges(start = as.integer(sj$start),
                                end   = as.integer(sj$end))
  )

  # ---- 2. GENCODE v29 gene models (lazy-loaded package data) -----------
  # utils::data() reads the installed package's lazy-load data db directly, so
  # this works whether or not prepTCGAdata has been library()'d, and — unlike
  # get(..., envir = asNamespace("prepTCGAdata")), which forces the package's
  # full namespace (and its heavy Imports: TCGAbiolinks, recount3, ...) to load
  # just to resolve one data object — it doesn't pull those in either. That's
  # the difference between this being instant and taking several seconds.
  genes_env <- new.env(parent = emptyenv())
  utils::data("gencode_v29_genes", package = "prepTCGAdata", envir = genes_env)
  genes_df <- genes_env$gencode_v29_genes

  # ---- 3. harmonise chromosome naming (chr1 vs 1) ----------------------
  sj_has_chr    <- any(grepl("^chr", as.character(sj[[chr_col]])))
  genes_has_chr <- any(grepl("^chr", genes_df$seqnames))
  if (sj_has_chr && !genes_has_chr) {
    genes_df$seqnames <- paste0("chr", genes_df$seqnames)
  } else if (!sj_has_chr && genes_has_chr) {
    genes_df$seqnames <- sub("^chr", "", genes_df$seqnames)
  }

  genes <- GenomicRanges::GRanges(
    seqnames = genes_df$seqnames,
    ranges   = IRanges::IRanges(start = genes_df$start, end = genes_df$end)
  )

  # ---- 4. overlap and collapse per interval --------------------------
  hits <- GenomicRanges::findOverlaps(sj_gr, genes, ignore.strand = TRUE)

  ov <- data.table::data.table(
    q         = S4Vectors::queryHits(hits),
    gene_id   = genes_df$gene_id[S4Vectors::subjectHits(hits)],
    gene_name = genes_df$gene_name[S4Vectors::subjectHits(hits)]
  )
  # Collapse per interval. De-duplicate on `gene_id` (GENCODE's key: each gene_id
  # has exactly one gene_name) and take the name from the *same* surviving rows,
  # so `gencode_gene_id` and `gencode_gene_name` always have the same number of
  # comma-separated entries, gene-for-gene in the same order. A gene_name shared
  # by several distinct loci (e.g. the Rfam-accession small RNAs like `RF00019`)
  # therefore appears once per overlapping gene_id rather than being collapsed to
  # one, which is what used to make the two lists disagree in length.
  agg <- ov[, {
    k <- !duplicated(gene_id)
    list(id   = paste(gene_id[k],   collapse = ","),
         name = paste(gene_name[k], collapse = ","))
  }, by = q]

  sj$gencode_gene_id   <- NA_character_
  sj$gencode_gene_name <- NA_character_
  sj$gencode_gene_id[agg$q]   <- agg$id
  sj$gencode_gene_name[agg$q] <- agg$name

  sj
}
