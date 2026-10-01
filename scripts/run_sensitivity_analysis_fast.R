run_sensitivity_analysis_fast <- function(
    exposomicset,
    base_formula,
    abundance_col = "counts",
    methods = c("limma_trend", "limma_voom", "DESeq2", "edgeR_quasi_likelihood"),
    scaling_methods = c("none", "TMM", "quantile"),
    contrasts = NULL,
    covariates_to_remove = NULL,
    da_res = NULL,
    pval_col = "adj.P.Val",
    logfc_col = "logFC",
    pval_threshold = 0.05,
    action = "add",
    bootstrap_n = 1
) {
  # ---- internal helpers ----
  
  .build_model_list <- function(base_formula, covariates_to_remove) {
    base_terms <- all.vars(base_formula)
    model_list <- list("full model" = base_formula)
    if (!is.null(covariates_to_remove)) {
      for (covar in covariates_to_remove) {
        reduced_terms <- setdiff(base_terms, covar)
        if (length(reduced_terms) > 1) {
          model_list[[paste("without", covar)]] <- as.formula(
            paste("~", paste(reduced_terms, collapse = " + "))
          )
        }
      }
    }
    model_list
  }
  
  # [Change #2] Resamples from pre-built exp_cache rather than the full MAE
  .resample_cache <- function(exp_cache) {
    purrr::map(exp_cache, function(se) {
      all_ids <- colnames(se)
      se[, sample(all_ids, replace = TRUE)]
    })
  }
  
  .filter_mae_to_da_res <- function(exposomicset, da_res) {
    sig_features <- da_res |> dplyr::distinct(feature, exp_name)
    for (en in unique(sig_features$exp_name)) {
      features_to_keep <- sig_features |>
        dplyr::filter(exp_name == en) |>
        dplyr::pull(feature)
      current_exp <- MultiAssayExperiment::experiments(exposomicset)[[en]]
      MultiAssayExperiment::experiments(exposomicset)[[en]] <- current_exp[
        rownames(current_exp) %in% features_to_keep,
      ]
    }
    message(sprintf(
      "Filtered MAE to %d features across %d experiments from supplied da_res.",
      nrow(sig_features), length(unique(sig_features$exp_name))
    ))
    exposomicset
  }
  
  .run_grid_cell_signs <- function(exp, exp_name, formula, method, scaling,
                                   abundance_col, contrasts, pval_col, logfc_col,
                                   pval_threshold) {
    if (method == "DESeq2") {
      SummarizedExperiment::assay(exp, abundance_col) <- round(
        SummarizedExperiment::assay(exp, abundance_col), 0
      )
    }
    
    res <- tryCatch(
      tidyexposomics:::.run_se_differential_abundance(
        se             = exp,
        formula        = formula,
        abundance_col  = abundance_col,
        method         = method,
        scaling_method = scaling,
        contrasts      = contrasts
      ),
      error = function(e) NULL
    )
    
    if (is.null(res)) return(NULL)
    
    res |>
      dplyr::select(
        feature = dplyr::any_of("feature"),
        logfc   = dplyr::all_of(logfc_col),
        pval    = dplyr::all_of(pval_col)
      ) |>
      dplyr::mutate(
        exp_name = exp_name,
        sign     = dplyr::case_when(
          is.na(.data$logfc) ~ NA_real_,
          .data$logfc > 0    ~  1,
          .data$logfc < 0    ~ -1,
          TRUE               ~  0
        ),
        is_sig = !is.na(.data$pval) & .data$pval < pval_threshold
      ) |>
      dplyr::select(feature, exp_name, sign, is_sig)
  }
  
  # [Change #1] Grid cells are now parallelised with future_pmap_dfr;
  # [Change #2] accepts a pre-built exp_cache instead of the full MAE
  .run_sign_grid <- function(exp_cache, model_list, methods, scalings,
                             abundance_col, contrasts, pval_col, logfc_col,
                             pval_threshold) {
    expand.grid(
      model_name = names(model_list),
      method     = methods,
      scaling    = scalings,
      exp_name   = names(exp_cache),
      stringsAsFactors = FALSE
    ) |>
      furrr::future_pmap_dfr(
        function(model_name, method, scaling, exp_name) {
          .run_grid_cell_signs(
            exp            = exp_cache[[exp_name]],
            exp_name       = exp_name,
            formula        = model_list[[model_name]],
            method         = method,
            scaling        = scaling,
            abundance_col  = abundance_col,
            contrasts      = contrasts,
            pval_col       = pval_col,
            logfc_col      = logfc_col,
            pval_threshold = pval_threshold
          )
        },
        # [Change #4] packages declared here so workers load them once at fork
        .options = furrr::furrr_options(
          seed     = TRUE,
          packages = c("MultiAssayExperiment", "SummarizedExperiment", "tidyexposomics")
        )
      ) |>
      dplyr::summarise(
        sign_consistency = ifelse(
          all(is.na(sign)), NA_real_,
          abs(sum(sign, na.rm = TRUE)) / sum(!is.na(sign))
        ),
        prop_sig = sum(is_sig, na.rm = TRUE) / dplyr::n(),
        .by = c(feature, exp_name)
      )
  }
  
  # ---- main logic ----
  
  model_list <- .build_model_list(base_formula, covariates_to_remove)
  
  if (!is.null(da_res)) {
    exposomicset <- .filter_mae_to_da_res(exposomicset, da_res)
  }
  
  # [Change #2] Build exp_cache once, outside the bootstrap loop
  exp_cache <- names(MultiAssayExperiment::experiments(exposomicset)) |>
    purrr::set_names() |>
    purrr::map(\(nm) tidyexposomics:::.update_assay_colData(exposomicset, nm))
  
  future::plan(future::multisession, workers = parallel::detectCores() - 1)
  
  if (bootstrap_n > 0) {
    feature_stability_df <- furrr::future_map(
      seq_len(bootstrap_n),
      # [Change #2] resample from cache; no library() calls needed here
      function(b) {
        .resample_cache(exp_cache) |>
          .run_sign_grid(
            model_list, methods, scaling_methods,
            abundance_col, contrasts, pval_col, logfc_col, pval_threshold
          ) |>
          dplyr::mutate(bootstrap_id = b)
      },
      # [Change #4] packages loaded on workers once at fork
      .options = furrr::furrr_options(
        seed     = TRUE,
        packages = c("MultiAssayExperiment", "SummarizedExperiment", "tidyexposomics")
      )
    ) |>
      purrr::list_rbind() |>
      dplyr::summarise(
        sign_consistency = mean(sign_consistency, na.rm = TRUE),
        prop_sig         = mean(prop_sig, na.rm = TRUE),
        .by = c(feature, exp_name)
      )
    
  } else {
    feature_stability_df <- .run_sign_grid(
      exp_cache, model_list, methods, scaling_methods,
      abundance_col, contrasts, pval_col, logfc_col, pval_threshold
    )
  }
  
  future::plan(future::sequential)
  
  if (action == "add") {
    all_metadata <- MultiAssayExperiment::metadata(exposomicset)
    all_metadata$differential_analysis$sensitivity_analysis <- list(
      feature_stability = feature_stability_df
    )
    MultiAssayExperiment::metadata(exposomicset) <- all_metadata
    
    MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(exposomicset)$summary$steps,
      list(run_sensitivity_analysis = list(
        timestamp = Sys.time(),
        params = list(
          formula              = format(base_formula),
          methods              = methods,
          scaling_methods      = scaling_methods,
          abundance_col        = abundance_col,
          contrasts            = contrasts,
          covariates_to_remove = covariates_to_remove,
          pval_threshold       = pval_threshold,
          bootstrap_n          = bootstrap_n,
          da_res_supplied      = !is.null(da_res)
        )
      ))
    )
    
    return(exposomicset)
    
  } else if (action == "get") {
    return(feature_stability_df)
  } else {
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
# Old code ----------
# run_sensitivity_analysis_fast <- function(
#     exposomicset,
#     base_formula,
#     abundance_col = "counts",
#     methods = c("limma_trend", "limma_voom", "DESeq2", "edgeR_quasi_likelihood"),
#     scaling_methods = c("none", "TMM", "quantile"),
#     contrasts = NULL,
#     covariates_to_remove = NULL,
#     pval_col = "adj.P.Val",
#     logfc_col = "logFC",
#     pval_threshold = 0.05,
#     logFC_threshold = log2(1),
#     score_thresh = NULL,
#     score_quantile = 0.9,
#     stability_metric = "stability_score",
#     action = "add",
#     bootstrap_n = 1
# ) {
#   # ---- internal helpers ----
#   
#   .build_model_list <- function(base_formula, covariates_to_remove) {
#     base_terms <- all.vars(base_formula)
#     model_list <- list("full model" = base_formula)
#     if (!is.null(covariates_to_remove)) {
#       for (covar in covariates_to_remove) {
#         reduced_terms <- setdiff(base_terms, covar)
#         if (length(reduced_terms) > 1) {
#           reduced_formula <- as.formula(
#             paste("~", paste(reduced_terms, collapse = " + "))
#           )
#           model_list[[paste("without", covar)]] <- reduced_formula
#         }
#       }
#     }
#     model_list
#   }
#   
#   .resample_MAE <- function(mae) {
#     all_ids <- colnames(mae) |> unlist() |> unique()
#     sample_ids <- sample(all_ids, replace = TRUE)
#     mae[, sample_ids]
#   }
#   
#   .run_sensitivity_grid <- function(exposomicset, model_list, methods, scalings, abundance_col, contrasts) {
#     exp_cache <- names(MultiAssayExperiment::experiments(exposomicset)) |>
#       purrr::set_names() |>
#       purrr::map(\(exp_name) tidyexposomics:::.update_assay_colData(exposomicset, exp_name))
#     
#     grid <- expand.grid(
#       model_name = names(model_list),
#       method     = methods,
#       scaling    = scalings,
#       exp_name   = names(MultiAssayExperiment::experiments(exposomicset)),
#       stringsAsFactors = FALSE
#     )
#     
#     purrr::pmap_dfr(grid, function(model_name, method, scaling, exp_name) {
#       exp <- exp_cache[[exp_name]]
#       
#       if (method == "DESeq2") {
#         SummarizedExperiment::assay(exp, abundance_col) <- round(
#           SummarizedExperiment::assay(exp, abundance_col), 0
#         )
#       }
#       
#       res <- .run_da_pipeline(exp, model_list[[model_name]], method, scaling, abundance_col, contrasts)
#       if (is.null(res)) return(NULL)
#       
#       res$model    <- model_name
#       res$exp_name <- exp_name
#       res
#     }) |>
#       tibble::as_tibble()
#   }
#   
#   .run_da_pipeline <- function(se, formula, method, scaling, abundance_col, contrasts) {
#     invisible(tidyexposomics:::.run_se_differential_abundance(
#       se            = se,
#       formula       = formula,
#       abundance_col = abundance_col,
#       method        = method,
#       scaling_method = scaling,
#       contrasts     = contrasts
#     ))
#   }
#   
#   .summarize_stable_features <- function(feature_stability_df, score_thresh, stability_metric) {
#     sum_df <- feature_stability_df |>
#       dplyr::count(
#         exp_name,
#         above = !!rlang::sym(stability_metric) > score_thresh
#       ) |>
#       tidyr::pivot_wider(
#         names_from  = above,
#         values_from = n,
#         values_fill = 0
#       ) |>
#       dplyr::rename(n_above = `TRUE`, n_below = `FALSE`) |>
#       dplyr::mutate(n = n_above + n_below)
#     
#     message("Number Of Features Above Threshold Of ", round(score_thresh, 2), ":")
#     message("----------------------------------------")
#     purrr::walk2(sum_df$exp_name, paste0(sum_df$n_above, "/", sum_df$n), message)
#   }
#   
#   # ---- main logic ----
#   
#   model_list <- .build_model_list(base_formula, covariates_to_remove)
#   
#   if (bootstrap_n > 0) {
#     future::plan(multisession, workers = parallel::detectCores() - 1)
#     
#     stability_list <- furrr::future_map(
#       seq_len(bootstrap_n),
#       function(b) {
#         message(sprintf("Running bootstrap iteration %d of %d", b, bootstrap_n))
#         .resample_MAE(exposomicset) |>
#           (\(mae_b) .run_sensitivity_grid(
#             mae_b, model_list, methods,
#             scaling_methods, abundance_col, contrasts
#           ))() |>
#           tidyexposomics:::.calculate_feature_stability(pval_col, logfc_col, pval_threshold) |>
#           dplyr::mutate(bootstrap_id = b)
#       },
#       .options = furrr::furrr_options(seed = TRUE)
#     )
#     
#     future::plan(sequential)
#     
#     feature_stability_df <- purrr::list_rbind(stability_list)
#     
#   } else {
#     feature_stability_df <- .run_sensitivity_grid(
#       exposomicset, model_list, methods,
#       scaling_methods, abundance_col, contrasts
#     ) |>
#       tidyexposomics:::.calculate_feature_stability(pval_col, logfc_col, pval_threshold)
#   }
#   
#   if (!stability_metric %in% colnames(feature_stability_df)) {
#     stop(sprintf("Invalid stability_metric: '%s'.", stability_metric))
#   }
#   
#   if (is.null(score_thresh)) {
#     score_thresh <- quantile(
#       feature_stability_df[[stability_metric]][feature_stability_df[[stability_metric]] > 0],
#       score_quantile,
#       na.rm = TRUE
#     )
#   }
#   
#   .summarize_stable_features(feature_stability_df, score_thresh, stability_metric)
#   
#   if (action == "add") {
#     all_metadata <- MultiAssayExperiment::metadata(exposomicset)
#     all_metadata$differential_analysis$sensitivity_analysis <- list(
#       feature_stability = feature_stability_df,
#       score_thresh      = score_thresh
#     )
#     MultiAssayExperiment::metadata(exposomicset) <- all_metadata
#     
#     step_record <- list(
#       run_sensitivity_analysis = list(
#         timestamp = Sys.time(),
#         params = list(
#           formula              = format(base_formula),
#           methods              = methods,
#           scaling_methods      = scaling_methods,
#           abundance_col        = abundance_col,
#           contrasts            = contrasts,
#           covariates_to_remove = covariates_to_remove,
#           pval_threshold       = pval_threshold,
#           logFC_threshold      = logFC_threshold,
#           score_quantile       = score_quantile,
#           stability_metric     = stability_metric,
#           bootstrap_n          = bootstrap_n
#         ),
#         notes = paste(
#           "Ran sensitivity analysis across", bootstrap_n,
#           "bootstrap iterations and", length(methods),
#           "methods and", length(scaling_methods), "scaling strategies.",
#           if (!is.null(covariates_to_remove)) {
#             paste("Covariates removed:", paste(covariates_to_remove, collapse = ", "))
#           } else {
#             "No covariates removed from base model."
#           }
#         )
#       )
#     )
#     
#     MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
#       MultiAssayExperiment::metadata(exposomicset)$summary$steps,
#       step_record
#     )
#     
#     return(exposomicset)
#     
#   } else if (action == "get") {
#     return(list(
#       feature_stability = feature_stability_df,
#       score_thresh      = score_thresh
#     ))
#   } else {
#     stop("Invalid action. Use 'add' or 'get'.")
#   }
# }

