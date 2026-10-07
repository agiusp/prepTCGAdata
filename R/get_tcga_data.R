#' get_tcga_data
#' @description Downloads TCGA splice junction data (GENCODE v29) and meta data,
#' annotates the junctions against the packaged GENCODE v29 gene models, scores
#' novel junctions, and prepares the sample meta data (survival and MSI status)
#' for overall-survival analysis.
#' @param ddir directory where you want the data to be written; the working
#' directory is changed to this location for the duration of the call and
#' restored on exit.
#' @param cohort TCGA cohort code to download (e.g. `"COAD"` or `"TCGA-COAD"`;
#' any `TCGA`, `-` and `_` characters are stripped).
#'
#' @details
#' Splice junctions and their per-sample counts are pulled from `recount3`
#' (`create_rse()`, `type = "jxn"`, `annotation = "gencode_v29"`) and written
#' out as `junction_counts.rds` (the raw junction x sample count matrix) and
#' `junction_metadata.rds` (one row per junction: coordinates, strand, and
#' the recount3/Snaptron annotation columns - `annotated`, `left_motif`,
#' `right_motif`, `left_annotated`, `right_annotated`).
#'
#' `junction_metadata.rds` is then overwritten with the result of
#' [annotate_sj()], which adds `gencode_gene_id` / `gencode_gene_name` columns
#' by overlapping every junction against the GENCODE release 29 gene models
#' shipped with the package ([gencode_v29_genes]) - release 29 because
#' recount3 builds its TCGA junctions against `"gencode_v29"`, so no GENCODE
#' file needs to be downloaded or supplied. That annotated table then drives
#' two further matrices, both restricted to novel (unannotated) junctions:
#' [count_novel_sjs()] tabulates novel-junction counts per gene and sample
#' (`novel_junction_counts_per_gene.rds`), and [novel_junction_RRS()] scores
#' each novel junction's relative read support against its gene's annotated
#' junctions (`novel_junction_RRS_scores.rds`). [novel_junction_ratios()] then
#' computes the fraction of expressed junctions that are novel, per gene
#' (`novel_junction_ratios.rds`) and per sample (added to `sample_metadata.rds`
#' as column `NJR`). The first three steps iterate over
#' the full junction set and can take a while (minutes to tens of minutes,
#' depending on cohort size).
#'
#' The novel-junction-per-gene counts are then rolled up by
#' [count_novel_sjs_per_pathway()] into scores for two gene-set collections:
#' every \pkg{xCell} cell-type signature (rows prefixed `"xCell:"`), and -
#' when the cohort has one - its \pkg{ConsensusTME} tumour-microenvironment
#' gene set (rows prefixed `"ConsTME:"`, e.g. `TCGA_COAD_...` uses
#' `ConsensusTME::consensusGeneSets[["COAD"]]`). The two are row-bound into a
#' single `novel_junction_counts_per_pathway.rds` matrix.
#'
#' `sample_metadata.rds` starts from recount3's own per-sample columns, then
#' has two further tables joined in, both keyed to the patient/sample via GDC
#' clinical and biospecimen records (`TCGAbiolinks::GDCquery_clinic()`):
#' overall-survival fields (`age_at_diagnosis_years`, `survival_days`,
#' `survival_years`, `censorship_flag`, and related derived columns), and, via
#' an internal `tcgamsi()` helper that queries the cBioPortal API,
#' microsatellite-instability status and scores (`msi_status`,
#' `msisensor_score`, `mantis_score`, `msi_call`, `msi_call_source`) - joined
#' on the sample's 15-character TCGA barcode (`msi.sample_id`), since
#' cBioPortal's sample id omits the portion/plate suffix letter GDC's own
#' `specimen_id` carries. `sample_metadata.rds` carries no Group /
#' SurviverGroup labelling of its own: that stratification (Histology / Stage
#' / age-band grouping, and the Good/Poor survivor label) happens
#' interactively in the SJSurv tab of the web app, where the age bands can be
#' tuned to the cohort at hand rather than fixed at download time.
#' @return Invisibly returns `NULL`. Called for its side effect of writing the
#' following files to `ddir` (each prefixed with `TCGA_<cohort>_`):
#' `junction_counts.rds`, `junction_metadata.rds`, `sample_metadata.rds`,
#' `novel_junction_counts_per_gene.rds`, `novel_junction_counts_per_pathway.rds`,
#' `novel_junction_RRS_scores.rds` and `novel_junction_ratios.rds`.
#' @seealso [annotate_sj()], [count_novel_sjs()], [count_novel_sjs_per_pathway()],
#'   [novel_junction_RRS()], [novel_junction_ratios()], [gencode_v29_genes]
#' @importFrom GSEABase geneIds setName
#' @export
get_tcga_data<-function(ddir,cohort){
  start_time <- Sys.time()

  if(length(grep('TCGA',cohort))>0)
    cohort<-gsub('-','',gsub('_','',gsub('TCGA','',cohort)))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(ddir)

  #!/usr/bin/env Rscript
  # ------------------------------------------------------------------
  # tcga_gencode_info.R
  #
  # Pulls:
  #   1. The list of TCGA project names/IDs from the GDC API
  #   2. The list of available GENCODE release versions from the
  #      GENCODE FTP site
  #
  # Requires: httr, jsonlite (for GDC), rvest (for GENCODE page parsing)
  # Install missing packages with:
  #   install.packages(c("httr", "jsonlite", "rvest"))
  # ------------------------------------------------------------------
  
  # --------------------------------------------------------------
  # 1. Get TCGA project names/IDs from the GDC API
  # --------------------------------------------------------------
  get_tcga_projects <- function() {
    base_url <- "https://api.gdc.cancer.gov/projects"
    
    filters <- list(
      op = "=",
      content = list(
        field = "program.name",
        value = "TCGA"
      )
    )
    
    params <- list(
      filters = toJSON(filters, auto_unbox = TRUE),
      fields = "project_id,name,primary_site,disease_type",
      format = "JSON",
      size = "200"
    )
    
    resp <- GET(base_url, query = params)
    stop_for_status(resp)
    
    parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
    hits <- parsed$data$hits
    
    if (is.null(hits) || nrow(hits) == 0) {
      warning("No TCGA projects returned from the GDC API.")
      return(data.frame())
    }
    
    # Keep it tidy and sorted
    hits <- hits[order(hits$project_id), c("project_id", "name", "primary_site", "disease_type")]
    rownames(hits) <- NULL
    hits
  }
  
  # --------------------------------------------------------------
  # 2. Get available GENCODE release versions
  #    (scrapes the human GENCODE release history page)
  # --------------------------------------------------------------
  get_gencode_versions <- function(species = c("human", "mouse")) {
    species <- match.arg(species)
    
    url <- if (species == "human") {
      "https://www.gencodegenes.org/human/releases.html"
    } else {
      "https://www.gencodegenes.org/mouse/releases.html"
    }
    
    page <- tryCatch(read_html(url), error = function(e) NULL)
    
    if (is.null(page)) {
      warning("Could not reach GENCODE releases page. Check network access.")
      return(character())
    }
    
    # The releases table has a column with release numbers, e.g. "Release 46"
    tbl <- page %>% html_element("table") %>% html_table(fill = TRUE)
    
    if (is.null(tbl) || nrow(tbl) == 0) {
      warning("Could not parse GENCODE release table; page layout may have changed.")
      return(character())
    }
    
    tbl
  }
  
  # --------------------------------------------------------------
  # Run
  # --------------------------------------------------------------
  if (sys.nframe() == 0) {
    cat("=== TCGA Projects (from GDC API) ===\n")
    tcga_projects <- get_tcga_projects()
    print(tcga_projects)
    
    cat("\n=== Available GENCODE Versions (human) ===\n")
    gencode_versions <- get_gencode_versions("human")
    print(gencode_versions)
    
    # Optionally save results to disk
    # write.csv(tcga_projects, "tcga_projects.csv", row.names = FALSE)
    # write.csv(gencode_versions, "gencode_versions.csv", row.names = FALSE)
  }
  
  # ==================================================================
  # Script: Extract All TCGA Cancer RNA Junction Matrices via recount3
  # ==================================================================
  
  # 1. Install and load required packages
  # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  # if (!requireNamespace("recount3", quietly = TRUE)) BiocManager::install("recount3")
  # if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")
  
  # library(recount3)
  # library(SummarizedExperiment)
  
  message(paste("--- Processing Cohort: TCGA-", cohort, " ---", sep=""))
  
  # 2. Query available projects in recount3 to target the right TCGA collection
  human_projects <- available_projects(organism = "human")
  project_info <- subset(human_projects, project == cohort & project_type == "data_sources")
  
  # 3. Pull the specific Splice Junction dataset (type = "junction")
  # Note: This will download files locally. Ensure you have ample disk space.
  message("Downloading junction data from recount3 repository...")
  print(project_info)
  rse_junction <- create_rse(
    project_info, 
    type = "jxn", 
    annotation = "gencode_v29"
  )
  
  # 4. Extract the Raw Splicing Junction Count Matrix
  # Rows = Unique genomic junction coordinates (chr:start-end:strand)
  # Columns = Patient/Sample TCGA barcodes
  junction_matrix <- assay(rse_junction, "counts")
  
  # 5. Extract structural junction metadata (genomic coordinates, hosting genes)
  junction_metadata <- as.data.frame(rowRanges(rse_junction))
  
  # 6. Extract patient clinical/sample metadata
  sample_metadata <- as.data.frame(colData(rse_junction))
  
  # 7. Print matrix dimensions to check dataset size
  message(paste("Successfully built matrix! Dimensions (Junctions x Samples):", 
                nrow(junction_matrix), "x", ncol(junction_matrix)))
  
  # 8. Export matrices locally as compressed files to protect system storage
  message("Saving data matrices to disk...")
  saveRDS(junction_matrix, file = paste0("TCGA_", cohort, "_junction_counts.rds"))
  saveRDS(junction_metadata, file = paste0("TCGA_", cohort, "_junction_metadata.rds"))
  saveRDS(sample_metadata, file = paste0("TCGA_", cohort, "_sample_metadata.rds"))
  
  #-----------------------------
  # 8.1. Annotate junctions by gene ids/names and count novel junctions per gene
  # against the packaged GENCODE v29 gene models. Both steps iterate over the
  # full junction set and can take a while.

  message("Annotating splice junctions ... this may take a few minutes ...")
  sjmeta<-annotate_sj(paste0("TCGA_", cohort, "_junction_metadata.rds"))
  saveRDS(sjmeta, file = paste0("TCGA_", cohort, "_junction_metadata.rds"))

  message("Counting novel splice junctions per gene ... this may take a while ...")
  cnsj<-count_novel_sjs(junction_matrix,sjmeta)
  saveRDS(cnsj, file = paste0("TCGA_", cohort, "_novel_junction_counts_per_gene.rds"))
  
  #collating count of novel sjs by pathway
  sigs <- xCell::xCell.data$signatures
  gsets<-lapply(sigs,function(x) geneIds(x))
  names(gsets) <- sapply(sigs, setName)
  cnsj.pway<-count_novel_sjs_per_pathway(cnsj,gsets)
  rownames(cnsj.pway)<-paste0('xCell:',rownames(cnsj.pway))
  if(length(intersect(cohort,names(ConsensusTME::consensusGeneSets)))>0){
    cnsj.pway2<-count_novel_sjs_per_pathway(cnsj,ConsensusTME::consensusGeneSets[[cohort]])
    rownames(cnsj.pway2)<-paste0('ConsTME:',rownames(cnsj.pway2))
    cnsj.pway<-rbind(cnsj.pway,cnsj.pway2)
  }
  saveRDS(cnsj.pway, file = paste0("TCGA_", cohort, "_novel_junction_counts_per_pathway.rds"))
  
  
  # 8.1. Compute RRS scores
  njrrs<-novel_junction_RRS (junction_matrix, sjmeta)
  saveRDS(njrrs,file = paste0("TCGA_", cohort, "_novel_junction_RRS_scores.rds"))

  # 8.2. Compute novel junction ratios: per gene (saved to file) and global
  # (joined onto sample_metadata as column NJR below)
  message("Computing novel junction ratios ...")
  njr.gene<-novel_junction_ratios(junction_matrix, sjmeta, global=FALSE)
  saveRDS(njr.gene,file = paste0("TCGA_", cohort, "_novel_junction_ratios.rds"))
  njr<-novel_junction_ratios(junction_matrix, sjmeta, global=TRUE)
  
  message("--- Process Finished! All RNA junction profiles exported successfully. ---")
  
  
  # =============================================================
  # Script: Get TCGA sample meta data and prepare for OS analysis
  # =============================================================
  
  # 1. Load meta data that has already been downloaded; we will join with OS
  #meta<-readRDS(lof[grep('sample_metadata.rds',lof)])
  meta<-sample_metadata
  meta$sample_id<-meta$tcga.gdc_file_id
  
  # 2. Query clinical and biospecimen data
  project_id<-paste0('TCGA-',cohort)
  clinical_raw <- GDCquery_clinic(project = project_id, type = "clinical")
  biospecimen_raw <- GDCquery_clinic(project = project_id, type = "biospecimen")
  
  # 3. Clean and map specimen metadata
  # In GDC biospecimen table, 'submitter_id' contains the sample barcode
  specimen_df <- biospecimen_raw %>%
    rename(specimen_id = submitter_id) %>%
    mutate(
      # Extract 12-character patient ID (e.g., TCGA-AA-3672)
      patient = substr(specimen_id, 1, 12),
      days_to_sample_procurement = suppressWarnings(as.numeric(days_to_sample_procurement))
    ) %>%
    # Keep primary tumor specimens (-01 barcode suffix)
    filter(grepl("-01[A-Z]?$", specimen_id)) %>%
    distinct(patient, .keep_all = TRUE) %>%
    select(patient, specimen_id, days_to_sample_procurement)
  
  # 4. Compute survival timeline and construct the output matrix
  survival_matrix <- clinical_raw %>%
    select(
      patient = submitter_id,
      age_at_diagnosis,
      vital_status,
      days_to_death,
      days_to_last_follow_up
    ) %>%
    left_join(specimen_df, by = "patient") %>%
    mutate(
      # Coerce numeric values safely
      age_at_diagnosis = suppressWarnings(as.numeric(age_at_diagnosis)),
      days_to_death = suppressWarnings(as.numeric(days_to_death)),
      days_to_last_follow_up = suppressWarnings(as.numeric(days_to_last_follow_up)),
      
      # Normalize age to years (TCGA reports age in days)
      age_at_diagnosis_years = ifelse(age_at_diagnosis > 120, age_at_diagnosis / 365.25, age_at_diagnosis),
      
      # 1 = Event (Death), 0 = Censored (Alive / Lost to follow-up)
      censorship_flag = ifelse(tolower(vital_status) == "dead", 1, 0),
      
      # Survival duration relative to diagnosis
      survival_days = ifelse(censorship_flag == 1, days_to_death, days_to_last_follow_up),
      survival_years = survival_days / 365.25,
      
      # Absolute patient age at event/censoring date
      age_at_event_or_censored = age_at_diagnosis_years + survival_years,
      
      # Patient age at specimen collection date
      procurement_offset_days = coalesce(days_to_sample_procurement, 0),
      specimen_collection_age = age_at_diagnosis_years + (procurement_offset_days / 365.25)
    ) %>%
    # Filter out missing baseline records
    filter(!is.na(survival_days), !is.na(age_at_diagnosis_years)) %>%
    select(
      patient_id = patient,
      specimen_id,
      age_at_diagnosis_years,
      survival_days,
      survival_years,
      age_at_event_or_censored,
      censorship_flag,
      specimen_collection_age
    ) %>%
    as.matrix()
  sm<-data.frame(survival_matrix)
  sm$tcga.cgc_case_id<-sm$patient_id
  
  #5. Join OS with meta and save. Cohort stratification (Histology / Stage /
  #   age band -> Group) and the Good/Poor SurviverGroup label are computed
  #   interactively in the SJSurv tab of the web app, not here.
  meta<-left_join(sm,meta)
  # saveRDS(meta,file = paste0("TCGA_", cohort, "_sample_metadata.rds"))
  
  # =============================================================
  # Script: Get TCGA MSI data
  # =============================================================
  
  tcgamsi <- function(COHORT_NAME,
                      derive = TRUE,
                      api = "https://www.cbioportal.org/api",
                      pancan_study = NULL,
                      pub_study = NULL,
                      quiet = FALSE) {
    
    if (!requireNamespace("httr", quietly = TRUE) ||
        !requireNamespace("jsonlite", quietly = TRUE)) {
      stop("tcgamsi() needs the packages 'httr' and 'jsonlite' - ",
           "install.packages(c('httr', 'jsonlite'))", call. = FALSE)
    }
    
    api <- sub("/+$", "", api)
    say <- function(...) if (!isTRUE(quiet)) message(sprintf(...))
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    # -- cBioPortal helpers ---------------------------------------------------- #
    
    # COAD and READ share cBioPortal's "coadread" study; everything else uses its
    # own project code lowercased.
    combined_slug <- function(code)
      switch(toupper(code), COAD = , READ = "coadread", tolower(code))
    
    api_get <- function(path, query = list()) {
      r <- httr::GET(paste0(api, path), query = query, httr::timeout(60),
                     httr::user_agent("tcgamsi() R function"))
      list(status = httr::status_code(r),
           body = tryCatch(
             jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                                simplifyDataFrame = TRUE),
             error = function(e) NULL))
    }
    
    study_exists <- function(study)
      api_get(paste0("/studies/", utils::URLencode(study, reserved = TRUE)))$status == 200
    
    # every SAMPLE- and PATIENT-level clinical value for a study, paged, as a
    # data.frame with sampleId / patientId / clinicalAttributeId / value.
    fetch_clinical <- function(study) {
      rows <- list()
      for (type in c("SAMPLE", "PATIENT")) {
        page <- 0L
        repeat {
          res <- api_get(
            sprintf("/studies/%s/clinical-data", utils::URLencode(study, reserved = TRUE)),
            list(clinicalDataType = type, projection = "SUMMARY",
                 pageSize = 100000L, pageNumber = page))
          if (res$status != 200 || is.null(res$body) || !length(res$body)) break
          df <- res$body
          if (is.null(df$sampleId)) df$sampleId <- NA_character_
          keep <- intersect(c("sampleId", "patientId", "clinicalAttributeId", "value"),
                            names(df))
          rows[[length(rows) + 1L]] <- df[, keep, drop = FALSE]
          if (nrow(df) < 100000L) break
          page <- page + 1L
        }
      }
      if (!length(rows)) return(NULL)
      do.call(rbind, rows)
    }
    
    # named vector (value keyed by sampleId) for the first attribute matching
    # `pattern`; a patient-level attribute is broadcast to that patient's samples.
    pick_attr <- function(clin, pattern) {
      hit <- unique(grep(pattern, clin$clinicalAttributeId, value = TRUE, ignore.case = TRUE))
      if (!length(hit)) return(NULL)
      sub <- clin[clin$clinicalAttributeId == hit[1L], , drop = FALSE]
      if (all(is.na(sub$sampleId) | sub$sampleId == "")) {
        smap <- unique(clin[!is.na(clin$sampleId) & clin$sampleId != "",
                            c("sampleId", "patientId")])
        sub <- merge(smap, sub[, c("patientId", "value")], by = "patientId")
      }
      stats::setNames(sub$value[!duplicated(sub$sampleId)],
                      sub$sampleId[!duplicated(sub$sampleId)])
    }
    
    norm_status <- function(x) {
      x <- toupper(trimws(x))
      out <- rep(NA_character_, length(x))
      out[x %in% c("MSI-H", "MSI_H", "MSIH", "HIGH")] <- "MSI-H"
      out[x %in% c("MSI-L", "MSI_L", "MSIL", "LOW")]  <- "MSI-L"
      out[x %in% c("MSS", "STABLE")]                  <- "MSS"
      out
    }
    
    # -- one cohort ---------------------------------------------------------- #
    
    one_cohort <- function(code) {
      code   <- toupper(trimws(code))
      pancan <- pancan_study %||% paste0(combined_slug(code), "_tcga_pan_can_atlas_2018")
      pub    <- pub_study    %||% paste0(combined_slug(code), "_tcga_pub")
      
      if (!study_exists(pancan))
        stop(sprintf("cBioPortal study '%s' not found - pass pancan_study= with the right id",
                     pancan), call. = FALSE)
      say("[%s] scores from %s", code, pancan)
      clin_pc <- fetch_clinical(pancan)
      if (is.null(clin_pc))
        stop(sprintf("no clinical data returned for '%s'", pancan), call. = FALSE)
      
      msisensor <- pick_attr(clin_pc, "^MSI_SENSOR_SCORE$")
      mantis    <- pick_attr(clin_pc, "^MSI_SCORE_MANTIS$")
      smap <- unique(clin_pc[!is.na(clin_pc$sampleId) & clin_pc$sampleId != "",
                             c("sampleId", "patientId")])
      
      msi_status <- NULL
      if (study_exists(pub)) {
        say("[%s] marker-panel MSI status from %s", code, pub)
        clin_pub <- fetch_clinical(pub)
        if (!is.null(clin_pub))
          msi_status <- pick_attr(clin_pub, "^MSI_STATUS(_[57]_MARKER_CALL)?$")
      } else {
        say("[%s] no %s study - msi_status left blank", code, pub)
      }
      
      g <- function(v, s) if (is.null(v)) rep(NA, length(s)) else unname(v[s])
      df <- data.frame(
        sample_id       = smap$sampleId,
        patient_id      = smap$patientId,
        cohort          = if (combined_slug(code) == "coadread") "COADREAD" else code,
        msi_status      = g(msi_status, smap$sampleId),
        msisensor_score = suppressWarnings(as.numeric(g(msisensor, smap$sampleId))),
        mantis_score    = suppressWarnings(as.numeric(g(mantis, smap$sampleId))),
        stringsAsFactors = FALSE)
      
      clean <- norm_status(df$msi_status)
      df$msi_status      <- ifelse(is.na(clean), df$msi_status, clean)  # keep 'Not Evaluable' etc.
      df$msi_call        <- clean
      df$msi_call_source <- ifelse(!is.na(clean), "marker_panel", NA_character_)
      if (isTRUE(derive)) {
        fill <- is.na(df$msi_call) & !is.na(df$mantis_score)
        df$msi_call[fill]        <- ifelse(df$mantis_score[fill] >= 0.4, "MSI-H", "MSS")
        df$msi_call_source[fill] <- "mantis>=0.4"
      }
      df[order(df$sample_id), , drop = FALSE]
    }
    
    # -- drive it ---------------------------------------------------------- #
    
    codes <- unlist(strsplit(paste(COHORT_NAME, collapse = ","), "[,[:space:]]+"))
    codes <- codes[nzchar(codes)]
    if (!length(codes)) stop("COHORT_NAME is empty", call. = FALSE)
    
    res <- do.call(rbind, lapply(codes, one_cohort))
    res <- res[!duplicated(res$sample_id), , drop = FALSE]
    rownames(res) <- NULL
    
    if (!isTRUE(quiet)) {
      tab <- table(factor(res$msi_call, c("MSI-H", "MSI-L", "MSS")), useNA = "ifany")
      message(sprintf("%d samples across %s", nrow(res),
                      paste(sort(unique(res$cohort)), collapse = ", ")))
      message("  msi_call: ",
              paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = "  "))
      message(sprintf("  MSIsensor: %d   MANTIS: %d   marker-panel status: %d",
                      sum(!is.na(res$msisensor_score)), sum(!is.na(res$mantis_score)),
                      sum(res$msi_call_source == "marker_panel", na.rm = TRUE)))
    }
    res
  }
  res<-tcgamsi(cohort) %>% dplyr::rename(msi.sample_id=sample_id)
  # cBioPortal's MSI sample id omits the portion/plate suffix letter GDC's
  # specimen_id carries (e.g. "TCGA-B0-5693-01" vs "TCGA-B0-5693-01A") -
  # truncate to the 15-character sample barcode before joining, or every
  # row silently fails to match and every msi.* column comes back NA.
  meta$msi.sample_id<-substr(meta$specimen_id,1,15)
  meta<-left_join(meta,res,by='msi.sample_id')
  
  #adding early/late Stage column
  q=grep('diagnoses.tumor_stage',colnames(meta),ignore.case = T)
  meta$Stage<-meta[,q]
  k<-union(grep('stage iii',meta$Stage),grep('stage iv',meta$Stage))
  meta$Stage[k]<-'Late Stage (3-4)'
  kk<-k
  k<-union(grep('stage i',meta$Stage),grep('stage ii',meta$Stage))
  meta$Stage[k]<-'Early Stage (1-4)'
  kk<-union(kk,k)
  meta$Stage[-kk]<-NA
  
  # meta$sample_id (set above from tcga.gdc_file_id) is a GDC file UUID, not
  # a sample barcode, so it won't generally match njrrs's (junction_matrix's)
  # column names - if it doesn't, look for a meta column whose values do
  # cover every njrrs sample barcode and use that instead.
  if(length(intersect(meta$sample_id,colnames(njrrs)))<ncol(njrrs)){
    n<-unlist(lapply(1:ncol(meta),
                     function(x) length(intersect(meta[,x],colnames(njrrs)))))
    q<-which(n==ncol(njrrs))[1]
    if(is.na(q)){
      stop("get_tcga_data(", cohort, "): no meta column's values cover all ",
           ncol(njrrs), " sample barcode(s) in njrrs/junction_matrix - ",
           "cannot determine sample_id.", call. = FALSE)
    }
    meta$sample_id<-meta[,q]
  }
  meta$NJR<-njr$NJR[match(meta$sample_id,njr$sample_id)]
  saveRDS(meta,file = paste0("TCGA_", cohort, "_sample_metadata.rds"))
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf("\nget_tcga_data(%s): completed %d sample(s) in %.2f seconds.\n",
              cohort, ncol(junction_matrix), elapsed))

}

