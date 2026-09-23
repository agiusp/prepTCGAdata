#' prepTCGAdata: Prepare TCGA Data for the SJV and SJVC Apps
#'
#' Downloads TCGA splice-junction data and clinical metadata via recount3 and
#' the GDC/cBioPortal APIs, annotates junctions against a GENCODE reference, and
#' prepares the metadata for overall-survival analysis.
#'
#' @keywords internal
#' @importFrom data.table := fread data.table
#' @importFrom dplyr %>% group_by summarize n full_join select rename mutate
#'   filter distinct left_join coalesce
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom httr GET stop_for_status content timeout user_agent status_code
#' @importFrom IRanges IRanges
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom recount3 available_projects create_rse
#' @importFrom rvest read_html html_element html_table
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay rowRanges colData
#' @importFrom TCGAbiolinks GDCquery_clinic
#' @importFrom utils txtProgressBar setTxtProgressBar URLencode
"_PACKAGE"

# Column names referenced through non-standard evaluation inside data.table []
# and dplyr verbs. Declaring them keeps R CMD check from reporting "no visible
# binding for global variable".
utils::globalVariables(c(
  ".",
  # data.table columns
  "type", "attribute", "gene_id", "gene_name", "seqnames", "start", "end",
  # subset() / dplyr columns
  "project", "project_type",
  "gencode_gene_id", "gencode_gene_name", "count", "count2",
  "submitter_id", "days_to_sample_procurement", "specimen_id", "patient",
  "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up",
  "age_at_diagnosis_years", "censorship_flag", "survival_days", "survival_years",
  "age_at_event_or_censored", "procurement_offset_days", "specimen_collection_age",
  "patient_id",
  # count_novel_sjs() / novel_junction_RRS() dplyr columns
  "annotated", "id", "name", "sj"
))
