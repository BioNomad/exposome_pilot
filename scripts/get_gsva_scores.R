get_gsva_scores <- function(
    se,
    gene_sets,
    prep_se   = FALSE,
    gene_symbol_col = "gene_symbol",
    log_trans = FALSE,
    assay     = "counts"
) {
  library(GSVA)
  
  if (prep_se) {
    se <- se |> (\(se) {
      assay(se, "counts") <- as.matrix(assay(se, "counts"))
      rowData(se) <- data.frame(gene_name = rownames(se))
      se
    })()
  }
  
  keep <- rownames(na.omit(assay(se, assay)))
  se   <- se[rownames(se) %in% keep, ]
  
  score_assay <- if (log_trans) {
    assay(se, "log") <- log2(assay(se, assay) + 1)
    "log"
  } else {
    assay
  }
  
  se <- se[rowData(se)[[gene_symbol_col]] %in% unique(unlist(gene_sets)), ]
  
  mat <- se |>
    pivot_se() |>
    dplyr::select(all_of(c(".sample", gene_symbol_col, score_assay))) |>
    dplyr::rename(score_val = all_of(score_assay)) |>
    filter(!is.na(score_val)) |>
    pivot_wider(names_from = gene_symbol_col,
                values_from = score_val,
                values_fn = mean) |>
    column_to_rownames(".sample") |>
    as.matrix()
  
  gsvaParam(t(mat), gene_sets, kcdf = "Gaussian", absRanking = TRUE) |>
    gsva() |>
    t() |>
    as.data.frame() |>
    rownames_to_column(".sample") |>
    left_join(se |> tidybulk::pivot_sample(), by = ".sample")
}