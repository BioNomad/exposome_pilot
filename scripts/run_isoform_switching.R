run_isoform_switching <- function(
    se, 
    assay="counts", 
    condition_col,
    covariates=NULL,
    condition_levels = NULL,
    gtf_path,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    verbose = TRUE,
    ignoreAfterPeriod    = TRUE,
    removeNonConvensionalChr = TRUE
) {
  # Load required packages
  suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
    library(rtracklayer)
    library(GenomicFeatures)
    library(tidyverse)
    library(edgeR)
    library(limma)
  })
  
  counts <- assay(se, assay)
  
  #---------------------------------
  #  Prepare counts + metadata
  #---------------------------------
  meta <- colData(se) %>% as.data.frame()
  
  # Filter samples with non-missing condition
  meta <- meta %>%
    filter(!is.na(.data[[condition_col]])) 
  
  
  
  if (!is.null(condition_levels)) {
    meta[[condition_col]] <- factor(meta[[condition_col]], levels = condition_levels)
  }
  
  meta$condition <- meta[[condition_col]]
  meta$sampleID <- rownames(meta)
  
  if(is.null(covariates)){
    meta <- meta |> 
      dplyr::select(sampleID,condition)
  }else{
    meta <- meta |> 
      dplyr::select(all_of(c("sampleID","condition",covariates)))
  }
  
  counts <- counts[, rownames(meta)]
  
  #---------------------------------
  #  Import into IsoformSwitchAnalyzeR
  #---------------------------------
  if (verbose) message("Creating IsoformSwitchAnalyzeR SwitchList...")
  
  switch_list <- importRdata(
    isoformCountMatrix   = counts,
    designMatrix         = meta,
    isoformExonAnnoation = gtf_path,
    ignoreAfterPeriod    = ignoreAfterPeriod,
    removeNonConvensionalChr = removeNonConvensionalChr,
    showProgress         = verbose
  )
  
  #---------------------------------
  #  Run Part 1 analysis
  #---------------------------------
  if (verbose) message("Running isoform switch analysis (Part 1)...")
  
  switch_part1 <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist   = switch_list,
    genomeObject         = genome,
    outputSequences      = FALSE,
    prepareForWebServers = FALSE
  )
  
  return(switch_part1)
}
