#' GENCODE v29 gene models
#'
#' The gene-level records of the GENCODE release 29 (GRCh38) annotation,
#' pre-parsed and shipped with the package so that [annotate_sj()],
#' [count_novel_sjs()] and [get_tcga_data()] can label splice junctions without
#' any GENCODE download or file path from the user.
#'
#' Release 29 is fixed on purpose: recount3 builds its TCGA splice-junction data
#' against `"gencode_v29"`, so the annotation reference must match.
#'
#' @format A data.frame with one row per gene and the columns:
#' \describe{
#'   \item{seqnames}{chromosome, GENCODE style (e.g. `"chr1"`).}
#'   \item{start, end}{1-based inclusive gene-body coordinates.}
#'   \item{strand}{`"+"` or `"-"`.}
#'   \item{gene_id}{versioned Ensembl/GENCODE gene id (e.g. `"ENSG00000223972.5"`).}
#'   \item{gene_name}{gene symbol (e.g. `"DDX11L1"`).}
#' }
#' The integer attribute `"gencode_release"` (29) records the source release.
#'
#' @source GENCODE release 29 GTF,
#'   \url{https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz}
#'   (see `data-raw/gencode_v29_genes.R`).
"gencode_v29_genes"
