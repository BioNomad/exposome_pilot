run_isoform_switching_satuRn <- function(
    se,
    counts_assay = "counts",
    tx2gene_gtf = NULL,               # optional: path to GTF for mapping
    tx2gene_df = NULL,                # OR provide your own data.frame(isoform_id, gene_id)
    formula = NULL,                   # model formula, e.g. ~ pftfev1fvc_actual + age + sex
    contrast_var = NULL,              # variable to test (string)
    min_count = 0,                    # filtering settings
    min_total_count = 0,
    parallel = FALSE,
    verbose = TRUE
) {
  # ---------------------------
  # 0. Basic input checks
  # ---------------------------
  stopifnot("SummarizedExperiment" %in% class(se))
  
  if (is.null(formula))
    stop("You must provide a model formula (e.g. ~ pft + age + sex).")
  
  if (is.null(contrast_var))
    stop("You must provide the variable name to test (contrast_var).")
  
  if (is.null(tx2gene_df) & is.null(tx2gene_gtf))
    stop("Provide either tx2gene_df or tx2gene_gtf to map transcripts to genes.")
  
  # ---------------------------
  # 1. Extract counts + metadata
  # ---------------------------
  counts <- SummarizedExperiment::assay(se, counts_assay)
  meta   <- as.data.frame(SummarizedExperiment::colData(se))
  
  # ---------------------------
  # 2. Build transcript-to-gene mapping
  # ---------------------------
  if (!is.null(tx2gene_gtf)) {
    if (verbose) message("Extracting tx→gene mapping from GTF...")
    
    suppressPackageStartupMessages({
      library(GenomicFeatures)
    })
    
    txdb <- txdbmaker::makeTxDbFromGFF(tx2gene_gtf)
    gtf_tx <- GenomicFeatures::transcripts(txdb, columns = c("transcript_id", "gene_id"))
    
    tx2gene <- as.data.frame(gtf_tx)[, c("transcript_id", "gene_id")]
    colnames(tx2gene) <- c("isoform_id", "gene_id")
    
  } else {
    tx2gene <- tx2gene_df
    colnames(tx2gene) <- c("isoform_id", "gene_id")
  }
  
  rownames(tx2gene) <- tx2gene$isoform_id
  
  # restrict to transcripts in counts
  tx2gene <- tx2gene[tx2gene$isoform_id %in% rownames(counts), ]
  
  # reorder to match count matrix
  tx2gene <- tx2gene[match(rownames(counts), tx2gene$isoform_id), ]
  
  # ---------------------------
  # 3. Remove single-isoform genes
  # ---------------------------
  suppressPackageStartupMessages({ library(dplyr) })
  
  tx2gene <- tx2gene %>%
    group_by(gene_id) %>%
    filter(n() > 1) %>%
    ungroup()
  
  counts <- counts[rownames(counts) %in% tx2gene$isoform_id, ]
  
  # reorder again
  tx2gene <- tx2gene[match(rownames(counts), tx2gene$isoform_id), ]
  
  if (verbose) message("Remaining transcripts: ", nrow(tx2gene), 
                       " across ", dplyr::n_distinct(tx2gene$gene_id), " genes.")
  
  # ---------------------------
  # 4. Filter low-expressed transcripts via edgeR
  # ---------------------------
  suppressPackageStartupMessages({ library(edgeR) })
  
  filt <- edgeR::filterByExpr(
    counts,
    min.count = min_count,
    min.total.count = min_total_count
  )
  
  counts <- counts[filt, ]
  tx2gene <- tx2gene[filt, ]
  
  # remove single-isoform genes again after filtering
  tx2gene <- tx2gene %>%
    group_by(gene_id) %>%
    filter(n() > 1) %>%
    ungroup()
  
  counts <- counts[rownames(counts) %in% tx2gene$isoform_id, ]
  
  tx2gene <- tx2gene[match(rownames(counts), tx2gene$isoform_id), ]
  
  if (verbose) message("After filtering: ", nrow(tx2gene), " transcripts remain.")
  
  # ---------------------------
  # 5. Build SummarizedExperiment for satuRn
  # ---------------------------
  suppressPackageStartupMessages({ library(SummarizedExperiment) })
  
  sumExp <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = meta,
    rowData = tx2gene
  )
  
  metadata(sumExp)$formula <- formula
  
  # ---------------------------
  # 6. Fit DTU models
  # ---------------------------
  suppressPackageStartupMessages({ library(satuRn) })
  
  if (verbose) message("Fitting satuRn DTU models...")
  
  sumExp <- satuRn::fitDTU(
    object = sumExp,
    formula = formula,
    parallel = parallel,
    verbose = verbose
  )
  
  # ---------------------------
  # 7. Build contrast matrix
  # ---------------------------
  design <- model.matrix(formula, data = meta)
  
  coef_names <- colnames(design)
  if (!(contrast_var %in% coef_names))
    stop("contrast_var '", contrast_var, "' not found in model coefficients.")
  
  L <- matrix(0, nrow = length(coef_names), ncol = 1,
              dimnames = list(coef_names, contrast_var))
  L[contrast_var, 1] <- 1
  
  # ---------------------------
  # 8. Test DTU
  # ---------------------------
  if (verbose) message("Testing DTU contrast: ", contrast_var)
  
  sumExp <- satuRn::testDTU(
    object = sumExp,
    contrasts = L,
    diagplot1 = TRUE,
    diagplot2 = TRUE,
    sort = TRUE
  )
  
  # ---------------------------
  # 9. Extract results table
  # ---------------------------
  result_name <- paste0("fitDTUResult_", contrast_var)
  res <- as.data.frame(SummarizedExperiment::rowData(sumExp)[[result_name]])
  res$isoform_id <- rownames(res)
  res <- res |> 
    left_join(rowData(se) |> 
                as.data.frame(),
              by=c("isoform_id"="transcript_id"))
  
  if (verbose) message("DTU completed successfully.")
  
  return(list(
    sumExp = sumExp,
    results = res,
    tx2gene = tx2gene,
    contrast = contrast_var
  ))
}
