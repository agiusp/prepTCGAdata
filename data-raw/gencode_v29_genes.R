# ---------------------------------------------------------------------------
# Build data/gencode_v29_genes.rda
#
# Pre-parses the GENCODE release 29 (GRCh38) annotation GTF down to the gene
# records that annotate_sj() / count_novel_sjs() need, and ships the result as
# package data so users never have to download or point at a GENCODE file.
#
# recount3 builds its TCGA splice junctions against "gencode_v29", so this
# release must not be changed without also changing the recount3 `annotation`
# argument in get_tcga_data().
#
# Source GTF (not redistributed with the package):
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
#
# Re-run with:  source("data-raw/gencode_v29_genes.R")
# ---------------------------------------------------------------------------

gtf_path <- path.expand("~/Work/SJ.Sep2026/Data/gencode.v29.annotation.gtf")
stopifnot(file.exists(gtf_path))

# Confirm the source really is release 29 before baking it in.
hdr <- readLines(gtf_path, n = 20L, warn = FALSE)
desc <- grep("^##description:", hdr, value = TRUE, ignore.case = TRUE)
ver  <- suppressWarnings(as.integer(sub(".*version[[:space:]]+([0-9]+).*", "\\1", desc[1])))
if (is.na(ver) || ver != 29L)
  stop("source GTF is not GENCODE release 29 (detected: ", ver, ")")

nskip <- sum(cumprod(startsWith(hdr, "#")))

gtf <- data.table::fread(
  gtf_path, sep = "\t", header = FALSE, quote = "", skip = nskip,
  showProgress = FALSE,
  col.names = c("seqnames", "source", "type", "start", "end",
                "score", "strand", "frame", "attribute")
)
gtf <- gtf[type == "gene"]
stopifnot(nrow(gtf) > 0L)

gencode_v29_genes <- data.frame(
  seqnames  = as.character(gtf$seqnames),
  start     = as.integer(gtf$start),
  end       = as.integer(gtf$end),
  strand    = as.character(gtf$strand),
  gene_id   = sub('.*gene_id "([^"]+)".*',   "\\1", gtf$attribute),
  gene_name = sub('.*gene_name "([^"]+)".*', "\\1", gtf$attribute),
  stringsAsFactors = FALSE
)

# GENCODE release the object represents; read back by the package for its
# provenance / sanity checks.
attr(gencode_v29_genes, "gencode_release") <- 29L

stopifnot(
  !anyNA(gencode_v29_genes$start),
  !anyNA(gencode_v29_genes$end),
  all(nzchar(gencode_v29_genes$gene_id)),
  all(nzchar(gencode_v29_genes$gene_name))
)

message(sprintf("gencode_v29_genes: %d genes, %d chromosomes",
                nrow(gencode_v29_genes),
                length(unique(gencode_v29_genes$seqnames))))

usethis::use_data(gencode_v29_genes, overwrite = TRUE, compress = "xz")
