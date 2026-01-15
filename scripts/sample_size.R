
calc_sample_sizes <- function(assoc_df, sd_y, alpha = 0.05, power = 0.90) {
  library(pwr)
  library(tidyverse)
  u <- 1  # testing the exposure
  
  assoc_df %>%
    mutate(
      beta = estimate,
      R2 = (beta^2) / (sd_y^2),
      f2  = R2 / (1 - R2),
      # If invalid, set NA
      f2  = ifelse(f2 <= 0 | is.infinite(f2), NA, f2),
      
      # Always return a list element
      pwr_result = map(f2, function(f2i) {
        if (is.na(f2i)) {
          return(NA)                 # list element = NA
        } else {
          return(pwr.f2.test(u = u, f2 = f2i, sig.level = alpha, power = power))
        }
      }),
      
      required_n = map_dbl(pwr_result, function(pr) {
        if (is.na(pr)[1]) return(NA_real_)  # safe test
        
        # extract v from pwr object and convert to total n
        return( ceiling(pr$v + u + 1) )
      })
    ) %>%
    select(term, beta, R2, f2, required_n)
}



calc_adj_sample_sizes <- function(assoc_df, alpha = 0.05, power = 0.90) {
  library(pwr)
  library(tidyverse)
  
  u <- 1  # testing one exposure
  
  assoc_df %>%
    mutate(
      t_value = statistic,
      df      = df, 
      R2_partial = (t_value^2) / (t_value^2 + df),
      f2 = R2_partial / (1 - R2_partial),
      
      pwr_result = map(f2, function(f2i) {
        if (is.na(f2i) || f2i <= 0) return(NA)
        pwr.f2.test(u = u, f2 = f2i, sig.level = alpha, power = power)
      }),
      
      required_n = map_dbl(pwr_result, function(pr) {
        if (is.na(pr)[1]) return(NA_real_)
        ceiling(pr$v + u + 1)
      })
    ) |> 
    select(term, estimate, statistic, df, R2_partial, f2, required_n)
}
