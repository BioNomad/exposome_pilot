library(dplyr)
library(purrr)
library(broom)
library(tibble)

# -------------------------------------------------------------
# Helpers
# -------------------------------------------------------------
.build_covariate_list <- function(covariates, covariates_to_remove = NULL) {
  covar_list <- list("full" = covariates)
  
  if (!is.null(covariates_to_remove)) {
    for (c in covariates_to_remove) {
      reduced <- setdiff(covariates, c)
      if (length(reduced) > 0)
        covar_list[[paste0("minus_", c)]] <- reduced
    }
  }
  
  covar_list
}

.bootstrap_resample <- function(df) {
  df[sample(seq_len(nrow(df)), replace = TRUE), , drop = FALSE]
}

.run_one_association <- function(df, exposure, outcome, covariates = character(), family = gaussian()) {
  bt  <- function(x) paste0("`", x, "`")
  rhs <- c(bt(exposure), bt(covariates)) |> paste(collapse = " + ")
  fml <- as.formula(paste(bt(outcome), "~", rhs))
  
  mod <- tryCatch(glm(fml, data = df, family = family), error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  
  out <- broom::tidy(mod) |> dplyr::filter(term == exposure)
  if (nrow(out) == 0) return(NULL)
  
  out |>
    dplyr::transmute(
      exposure  = exposure,
      estimate  = estimate,
      std.error = std.error,
      statistic = statistic,
      p_value   = p.value
    )
}

.calculate_stability <- function(sensitivity_df, pval_threshold = 0.05, logFC_threshold = 0) {
  sensitivity_df |>
    mutate(
      signif = (p_value <= pval_threshold) & (abs(estimate) >= logFC_threshold),
      sign   = sign(estimate)
    ) |>
    group_by(exposure) |>
    summarise(
      n_tests          = n(),
      n_signif         = sum(signif, na.rm = TRUE),
      prop_signif      = n_signif / n_tests,
      mean_effect      = mean(estimate, na.rm = TRUE),
      sd_effect        = sd(estimate, na.rm = TRUE),
      sign_consistency = ifelse(
        all(is.na(sign)), NA_real_,
        abs(sum(sign, na.rm = TRUE)) / n_tests
      ),
      stability_score  = prop_signif * sign_consistency,
      .groups = "drop"
    )
}

# -------------------------------------------------------------
# Main function
# -------------------------------------------------------------
run_exwas_sensitivity <- function(
    df,
    exposures,
    outcome,
    covariates,
    covariates_to_remove = NULL,
    family          = gaussian(),
    pval_threshold  = 0.05,
    logFC_threshold = 0,
    score_quantile  = 0.9,
    bootstrap_n     = 1,
    parallel        = FALSE
) {
  
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE))
      stop("Install the `furrr` package to use parallel = TRUE.")
    message("Running in parallel — ensure a future plan is set, e.g. future::plan(multisession).")
  }
  
  covar_list <- .build_covariate_list(covariates, covariates_to_remove)
  
  grid <- tidyr::expand_grid(
    bootstrap  = seq_len(bootstrap_n),
    model_name = names(covar_list),
    exposure   = exposures
  )
  
  worker <- function(i) {
    row  <- grid[i, ]
    covs <- covar_list[[row$model_name]]
    df_b <- .bootstrap_resample(df)
    
    res <- .run_one_association(df_b, row$exposure, outcome, covs, family)
    if (is.null(res)) return(NULL)
    
    res$model     <- row$model_name
    res$bootstrap <- row$bootstrap
    res
  }
  
  all_results <- if (parallel) {
    furrr::future_map_dfr(seq_len(nrow(grid)), worker, .progress = TRUE)
  } else {
    purrr::map_dfr(seq_len(nrow(grid)), worker)
  }
  
  feature_stability <- .calculate_stability(all_results, pval_threshold, logFC_threshold)
  
  score_thresh <- quantile(
    feature_stability$stability_score[feature_stability$stability_score > 0],
    score_quantile,
    na.rm = TRUE
  )
  
  list(
    sensitivity_df    = all_results,
    feature_stability = feature_stability,
    score_thresh      = score_thresh,
    stable_features   = feature_stability |> filter(stability_score >= score_thresh)
  )
}
