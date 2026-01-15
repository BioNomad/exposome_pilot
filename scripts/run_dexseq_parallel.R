run_dexseq_parallel <- function(
    se,
    condition_col,
    covariates = NULL,
    gtf_path,
    design_full = NULL,
    design_reduced = NULL,
    cores = 4,
    backend = c("Multicore", "Snow"),
    verbose = TRUE
) {
  suppressPackageStartupMessages({
    library(DEXSeq)
    library(GenomicFeatures)
    library(GenomicAlignments)
    library(SummarizedExperiment)
    library(BiocParallel)
    library(txdbmaker)
    library(dplyr)
    library(tibble)
  })
  
  backend <- match.arg(backend)
  
  #-----------------------------
  # PARALLEL BACKEND FIX
  #-----------------------------
  if (.Platform$OS.type == "windows" && backend == "Multicore") {
    message("⚠️ Windows does not support MulticoreParam. Switching to SnowParam automatically.")
    backend <- "Snow"
  }
  
  if (verbose) message("▶ Setting up parallel backend…")
  
  BPPARAM <- switch(
    backend,
    Multicore = MulticoreParam(workers = cores),
    Snow = SnowParam(workers = cores, type = "SOCK")
  )
  
  #-----------------------------
  # METADATA / DESIGN
  #-----------------------------
  if (verbose) message("▶ Extracting metadata…")
  meta <- as.data.frame(colData(se))
  
  meta$condition <- factor(meta[[condition_col]])
  
  if (!is.null(covariates)) {
    missing <- setdiff(covariates, names(meta))
    if (length(missing) > 0)
      stop("Missing covariates: ", paste(missing, collapse = ", "))
  }
  
  if (is.null(design_full))
    design_full <- ~ sample + exon + condition:exon
  
  if (is.null(design_reduced))
    design_reduced <- ~ sample + exon
  
  #-----------------------------
  # FLATTEN GTF (FIXED)
  #-----------------------------
  if (verbose) message("▶ Loading GTF and flattening annotation…")
  
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)
  
  flattened <- GenomicFeatures::exonicParts(
    txdb,
    linked.to.single.gene.only = TRUE
  )
  
  names(flattened) <- sprintf(
    "%s:E%03d", flattened$gene_id, flattened$exonic_part
  )
  
  #-----------------------------
  # MAP COUNTS → EXON BINS
  #-----------------------------
  if (verbose) message("▶ Summarizing counts into exon bins (parallel)…")
  
  if (!"gene_id" %in% colnames(rowData(se)))
    stop("rowData(se) must contain 'gene_id'.")
  
  gene_map <- rowData(se)$gene_id
  counts_mat <- assay(se)
  
  exon_gene_ids <- flattened$gene_id
  exon_ids <- names(flattened)
  
  exon_counts <- bplapply(
    seq_along(exon_ids),
    function(i) {
      g <- exon_gene_ids[i]
      idx <- which(gene_map == g)
      if (length(idx) > 0) {
        colSums(counts_mat[idx, , drop = FALSE])
      } else {
        rep(0, ncol(counts_mat))
      }
    },
    BPPARAM = BPPARAM
  )
  
  exon_count_matrix <- do.call(rbind, exon_counts)
  rownames(exon_count_matrix) <- exon_ids
  
  se_counts <- SummarizedExperiment(
    assays = list(counts = exon_count_matrix),
    rowRanges = flattened,
    colData = meta
  )
  
  #-----------------------------
  # DEXSEQ WORKFLOW
  #-----------------------------
  meta$sample <- factor(rownames(meta))
  
  if (verbose) message("▶ Constructing DEXSeqDataSetFromSE...")
  dxd <- DEXSeqDataSetFromSE(
    se_counts,
    design = design_full
  )
  
  if (verbose) message("▶ Estimating size factors…")
  dxd <- estimateSizeFactors(dxd)
  
  if (verbose) message("▶ Estimating dispersions (parallel)…")
  dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
  
  if (verbose) message("▶ Testing DEU (parallel)…")
  dxd <- testForDEU(
    dxd,
    reducedModel = design_reduced,
    fullModel = design_full,
    BPPARAM = BPPARAM
  )
  
  if (verbose) message("▶ Estimating exon fold changes (parallel)…")
  dxd <- estimateExonFoldChanges(
    dxd,
    fitExpToVar = "condition",
    BPPARAM = BPPARAM
  )
  
  if (verbose) message("▶ Producing results table…")
  dxr <- DEXSeqResults(dxd)
  
  return(list(
    dxd = dxd,
    results = dxr,
    flattened = flattened,
    parallel = BPPARAM
  ))
}
