pivot_se <- function(se, features = NULL, id_col = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  
  # Optional feature filtering
  if (!is.null(features)) {
    se <- se[rownames(se) %in% features, , drop = FALSE]
  }
  
  # Extract metadata
  rd <- as.data.frame(rowData(se))
  cd <- as.data.frame(colData(se))
  
  # Pivot each assay to long format
  assay_list <- lapply(names(assays(se)), function(aname) {
    a <- assays(se)[[aname]]
    tidyr::pivot_longer(
      tibble::as_tibble(a, rownames = ".feature"),
      cols = -.feature,
      names_to  = ".sample",
      values_to = aname
    )
  })
  
  # Merge assays
  assay_df <- Reduce(
    function(x, y) dplyr::left_join(x, y, by = c(".feature", ".sample")),
    assay_list
  )
  
  # Attach rowData and colData
  assay_df <- assay_df |>
    dplyr::left_join(
      tibble::rownames_to_column(rd, ".feature"),
      by = ".feature"
    ) |>
    dplyr::left_join(
      tibble::rownames_to_column(cd, ".sample"),
      by = ".sample"
    )
  
  # Optional identifier column (useful when binding multiple SEs)
  if (!is.null(id_col)) {
    assay_df <- dplyr::mutate(
      assay_df,
      se_id = id_col,
      .before = .feature
    )
  }
  
  assay_df
}
