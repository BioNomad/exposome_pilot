covar_screen <- function(
    data, 
    outcome,
    base_covars,
    candidate_covars,
    vif_cutoff = 8
){
  library(dplyr)
  library(broom)
  library(rsq)
  library(car)
  library(purrr)
  library(tibble)
  
  #-------------------------------------------
  # Safe VIF helper
  #-------------------------------------------
  safe_vif <- function(fit, candidate_covars) {
    out <- tibble(
      var = candidate_covars,
      VIF = NA_real_
    )
    
    if (is.null(fit)) return(out)
    
    vv <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (is.null(vv)) return(out)
    
    if (is.null(names(vv))) {
      names(vv) <- candidate_covars[seq_along(vv)]
    }
    
    tibble(
      var = names(vv),
      VIF = as.numeric(vv)
    ) %>% 
      filter(var %in% candidate_covars)
  }
  
  # --- Univariate tests ---------------
  # does the covariate associate with the outcome on its own?
  test_univariate <- function(var){
    f <- as.formula(paste(outcome, "~", var))
    fit <- lm(f, data = data)
    tibble(
      var = var,
      uni_p = glance(fit)$p.value,
      uni_adj_r2 = glance(fit)$adj.r.squared
    )
  }
  
  uni_results <- map_df(candidate_covars, test_univariate)
  
  # --- Incremental tests ----
  # does adding this covariate improve the prediction beyond the base
  # covariates
  base_formula <- as.formula(
    paste(outcome, "~", paste(base_covars, collapse = "+"))
  )
  base_fit <- lm(base_formula, data = data)
  
  test_incremental <- function(var){
    new_formula <- as.formula(
      paste(outcome, "~", paste(c(base_covars, var), collapse = "+"))
    )
    new_fit <- lm(new_formula, data = data)
    
    tibble(
      var = var,
      delta_adj_r2 = glance(new_fit)$adj.r.squared - glance(base_fit)$adj.r.squared,
      delta_AIC = AIC(new_fit) - AIC(base_fit),
      delta_BIC = BIC(new_fit) - BIC(base_fit)
    )
  }
  
  inc_results <- map_df(candidate_covars, test_incremental)
  
  # --- VIF with safe handler ----
  # Is this covariate redundant given the VIF?
  full_formula <- as.formula(
    paste(outcome, "~", paste(c(base_covars, candidate_covars), collapse="+"))
  )
  
  full_fit <- tryCatch(
    lm(full_formula, data = data),
    error = function(e) NULL
  )
  
  vif_df <- safe_vif(full_fit, candidate_covars)
  
  # --- Combine results -----------
  # does the covariate have a positivie change in adj r2
  # and is the vif below the cutoff, and is the change in 
  # AIC/BIC below 0 (so improving the fit)
  combined <- uni_results %>%
    left_join(inc_results, by = "var") %>%
    left_join(vif_df, by = "var") %>%
    arrange(desc(delta_adj_r2)) %>%
    mutate(
      recommended = case_when(
        delta_adj_r2 > 0 &
          delta_AIC < 0 &
          delta_BIC < 0 &
          (is.na(VIF) | VIF < vif_cutoff) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  list(
    univariate = uni_results,
    incremental = inc_results,
    vif = vif_df,
    combined_ranked = combined,
    recommended_covariates = combined %>% filter(recommended) %>% pull(var)
  )
}
