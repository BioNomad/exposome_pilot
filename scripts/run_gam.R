run_gam <- function(
    df,
    codebook=NULL,
    outcome,
    covariates,
    feature_set,
    family = gaussian(),
    k = 4
) {
  require(mgcv)
  res <- map(feature_set, ~ {
    complete_df <- df |>
      dplyr::select(all_of(c(outcome, covariates, .x))) |>
      tidyr::drop_na()
    
    if (length(unique(complete_df[[.x]])) < k) {
      message(sprintf("Skipping %s — fewer than %d unique values after dropping NAs", .x, k))
      return(NULL)
    }
    
    fm <- as.formula(paste(
      outcome, "~",
      paste0("s(", .x, ")"), "+",
      paste(covariates, collapse = " + ")
    ))
    
    fit <- tryCatch(
      gam(fm, data = complete_df, family = family),
      error = function(e) {
        message(sprintf("GAM failed for %s: %s", .x, e$message))
        NULL
      }
    )
    
    if (is.null(fit)) return(NULL)
    
    broom::tidy(fit) |>
      mutate(variable = .x)
    
  }) |>
    purrr::compact() |>
    bind_rows()
  
  if(is.null(codebook)){
    res
  } else{
    res <- res |> 
      left_join(
        codebook,
        by="variable"
      )
  }
}


