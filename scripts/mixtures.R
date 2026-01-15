run_bws_group <- function(group_name,
                          dat,
                          outcome,
                          covariates,
                          family = "gaussian") {
  
  library(bws)
  
  mix_vars <- group_list[[group_name]]
  
  X <- dat %>%
    select(all_of(mix_vars)) %>%
    as.matrix()
  
  Y <- dat[[outcome]]
  
  Z <- dat |> 
    select(all_of(covariates)) %>%
    mutate(across(where(is.factor), as.numeric)) %>%
    as.matrix()
  
  fit <- bws(
    iter = 3000,
    y = Y,
    X = X,
    Z = Z,
    family = family,
    chains = 4,
    cores = 4,
    refresh = 0
  )
  
  list(group = group_name, fit = fit, vars = mix_vars)
}

