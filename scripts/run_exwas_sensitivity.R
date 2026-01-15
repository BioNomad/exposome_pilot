library(dplyr)
library(purrr)
library(broom)
library(tibble)

#------------------------------------------
# Helper: Build list of covariate sets
#------------------------------------------
build_covariate_list <- function(covariates, covariates_to_remove = NULL) {
  covar_list <- list("full" = covariates)
  
  if (!is.null(covariates_to_remove)) {
    for (c in covariates_to_remove) {
      reduced <- setdiff(covariates, c)
      if (length(reduced) > 0) {
        covar_list[[paste0("minus_", c)]] <- reduced
      }
    }
  }
  
  covar_list
}

#------------------------------------------
# Helper: Resample rows (bootstrap)
#------------------------------------------
bootstrap_resample <- function(df) {
  df[sample(seq_len(nrow(df)), replace = TRUE), , drop = FALSE]
}

#------------------------------------------
# Run GLM association for one exposure & one covariate set
#------------------------------------------

run_one_association <- function(df,
                                exposure,
                                outcome,
                                covariates = character(),
                                family = gaussian()) {
  
  # helper: backtick variable names safely
  bt <- function(x) paste0("`", x, "`")
  
  rhs <- c(bt(exposure), bt(covariates)) |>
    paste(collapse = " + ")
  
  fml <- as.formula(paste(
    bt(outcome),
    "~",
    rhs
  ))
  
  mod <- tryCatch(
    glm(fml, data = df, family = family),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)
  
  tidy_res <- broom::tidy(mod)
  
  out <- tidy_res |>
    dplyr::filter(term == exposure)
  
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

# run_one_association <- function(df, exposure, outcome, covariates, family) {
#   
#   rhs <- c(exposure, covariates) |> paste(collapse = " + ")
#   formula <- as.formula(paste(outcome, "~", rhs))
#   
#   mod <- tryCatch(glm(formula, data = df, family = family), error = function(e) NULL)
#   if (is.null(mod)) return(NULL)
#   
#   tidy_res <- broom::tidy(mod)
#   out <- tidy_res |> filter(term == exposure)
#   
#   if (nrow(out) == 0) return(NULL)
#   
#   out$estimate <- out$estimate
#   out$p_value <- out$p.value
#   out$exposure <- exposure
#   out
# }

#------------------------------------------
# Calculate stability across all model runs
#------------------------------------------
calculate_stability <- function(sensitivity_df,
                                pval_threshold = 0.05,
                                logFC_threshold = 0) {
  
  sensitivity_df %>%
    mutate(
      signif = (p_value <= pval_threshold) & (abs(estimate) >= logFC_threshold),
      sign = sign(estimate)
    ) %>%
    group_by(exposure) %>%
    summarise(
      n_tests = n(),
      n_signif = sum(signif, na.rm = TRUE),
      prop_signif = n_signif / n_tests,
      mean_effect = mean(estimate, na.rm = TRUE),
      sd_effect = sd(estimate, na.rm = TRUE),
      sign_consistency =
        ifelse(all(is.na(sign)), NA_real_,
               abs(sum(sign, na.rm = TRUE)) / n_tests),
      stability_score = prop_signif * sign_consistency,
      .groups = "drop"
    )
}

#------------------------------------------
# Main LOCAL Sensitivity Function
#------------------------------------------
run_exwas_sensitivity <- function(
    df,
    exposures,
    outcome,
    covariates,
    covariates_to_remove = NULL,
    family = gaussian(),
    pval_threshold = 0.05,
    logFC_threshold = 0,
    score_quantile = 0.9,
    bootstrap_n = 1
) {
  
  covar_list <- build_covariate_list(covariates, covariates_to_remove)
  
  all_results <- map_dfr(seq_len(bootstrap_n), function(b) {
    
    message(sprintf("Bootstrap %d of %d", b, bootstrap_n))
    df_b <- bootstrap_resample(df)
    
    # Grid across covariate sets and exposures
    map_dfr(names(covar_list), function(model_name) {
      covs <- covar_list[[model_name]]
      
      map_dfr(exposures, function(exp){
        res <- run_one_association(df_b, exp, outcome, covs, family)
        if (is.null(res)) return(NULL)
        res$model <- model_name
        res$bootstrap <- b
        res
      })
    })
  })
  
  # Compute stability metrics
  feature_stability <- calculate_stability(
    all_results,
    pval_threshold = pval_threshold,
    logFC_threshold = logFC_threshold
  )
  
  # Determine threshold for stable features
  score_thresh <- quantile(
    feature_stability$stability_score[feature_stability$stability_score > 0],
    score_quantile,
    na.rm = TRUE
  )
  
  list(
    sensitivity_df = all_results,
    feature_stability = feature_stability,
    score_thresh = score_thresh,
    stable_features = feature_stability |>
      filter(stability_score >= score_thresh)
  )
}
