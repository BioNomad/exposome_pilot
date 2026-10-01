run_sva_limma_trend <- function(mae, assay_name, base_formula) {
  message("Processing: ", assay_name)
  
  se  <- tidyexposomics:::.update_assay_colData(mae, assay_name)
  mat <- SummarizedExperiment::assay(se, "counts") |> as.matrix()
  df  <- SummarizedExperiment::colData(se) |> as.data.frame()
  
  # --- SVA ---
  mod  <- model.matrix(base_formula, data = df)
  mod0 <- model.matrix(
    update(base_formula, ~ . - pftfev1fvc_actual), 
    data = df
  )
  
  svobj <- sva::sva(mat, mod, mod0)
  message("  n.sv = ", svobj$n.sv)
  
  # add SVs to colData
  sv_df <- svobj$sv |>
    as.data.frame() |>
    setNames(paste0("SV", seq_len(svobj$n.sv)))
  
  df_sv <- dplyr::bind_cols(df, sv_df)
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(df_sv)
  
  # build formula dynamically from detected SVs
  sv_terms <- paste0("SV", seq_len(svobj$n.sv)) |> paste(collapse = " + ")
  full_formula <- update(base_formula, paste("~ . +", sv_terms))
  message("  formula: ", deparse(full_formula))
  
  # --- limma-trend ---
  tidyexposomics:::.run_limma_trend(
    se             = se,
    formula        = full_formula,
    abundance_col  = "counts",
    scaling_method = "none"
  ) |>
    dplyr::mutate(exp_name = assay_name, n_sv = svobj$n.sv)
}

