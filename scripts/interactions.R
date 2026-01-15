# test_pair <- function(df, x, y) {
#   
#   # Additive model
#   f_add <- reformulate(
#     c(x, y, "gli_age", "gli_sex", "gli_height", "fis"),
#     response = "pftfev1fvc_actual"
#   )
#   
#   # Full interaction model
#   f_int <- reformulate(
#     c(paste0(x, "*", y), "gli_age", "gli_sex", "gli_height", "fis"),
#     response = "pftfev1fvc_actual"
#   )
#   
#   m_add <- glm(f_add, data = df)
#   m_int <- glm(f_int, data = df)
#   
#   # Extract interaction estimate + p-value
#   int_term <- broom::tidy(m_int) |> 
#     filter(term == paste0(x, ":", y))
#   
#   # LRT between nested models
#   LRT <- anova(m_add, m_int, test = "LRT")
#   LRT_p <- LRT$`Pr(>Chi)`[2]
#   
#   tibble(
#     var1 = x,
#     var2 = y,
#     int_est = int_term$estimate,
#     p_int = int_term$p.value,
#     AIC_add = AIC(m_add),
#     AIC_int = AIC(m_int),
#     AIC_supports_add = AIC(m_add) < AIC(m_int),
#     LRT_p = LRT_p,
#     
#     # Classify mode of joint action
#     mode = dplyr::case_when(
#       int_est < 0 & p_int < 0.05 ~ "Synergistic",
#       int_est > 0 & p_int < 0.05 ~ "Antagonistic",
#       (p_int >= 0.05) & (AIC(m_add) < AIC(m_int)) & (LRT_p > 0.05) ~ "Additive",
#       TRUE ~ "Other"
#     )
#   )
# }

# --- Updated Function ---------------
library(broom)
library(dplyr)
library(purrr)

test_pair <- function(df, x, y) {
  
  # ----- Model formulas -----
  f_add <- reformulate(
    c(x, y, "gli_age", "gli_sex", "gli_height", "fis"),
    response = "pftfev1fvc_actual"
  )
  
  f_int <- reformulate(
    c(paste0(x, "*", y), "gli_age", "gli_sex", "gli_height", "fis"),
    response = "pftfev1fvc_actual"
  )
  
  # ----- Fit models -----
  m_add <- glm(f_add, data = df)
  m_int <- glm(f_int, data = df)
  
  # ----- Extract coefficients -----
  coefs <- broom::tidy(m_int)
  
  beta1 <- coefs |> filter(term == x) |> pull(estimate)
  beta2 <- coefs |> filter(term == y) |> pull(estimate)
  
  int_term <- coefs |> filter(term == paste0(x, ":", y))
  
  int_est <- int_term$estimate
  p_int <- int_term$p.value
  
  # ----- Likelihood ratio test -----
  LRT <- anova(m_add, m_int, test = "LRT")
  LRT_p <- LRT$`Pr(>Chi)`[2]
  
  AIC_add <- AIC(m_add)
  AIC_int <- AIC(m_int)
  AIC_supports_add <- AIC_add < AIC_int
  
  #  Mode (Synergy / Antagonism / Additive)
  mode <- case_when(
    int_est < 0 & p_int < 0.05 ~ "Synergistic",
    int_est > 0 & p_int < 0.05 ~ "Antagonistic",
    p_int >= 0.05 & AIC_supports_add & LRT_p > 0.05 ~ "Additive",
    TRUE ~ "Other"
  )
  
  #  Interaction Direction (Positive or Negative)
  interaction_direction <- case_when(
    int_est < 0 ~ "Negative interaction",
    int_est > 0 ~ "Positive interaction",
    TRUE ~ "No direction"
  )
  
  #  Effect Directionality for FEV1/FVC (Higher = Better)
  # logic combining effects
  effect_direction <- case_when(
    
    # Synergy/Antagonism: both harmful exposures
    beta1 < 0 & beta2 < 0 & int_est < 0 ~ "Amplifies Negative",
    beta1 < 0 & beta2 < 0 & int_est > 0 ~ "Buffers Negative",
    
    # Synergy/Antagonism: both beneficial exposures
    beta1 > 0 & beta2 > 0 & int_est < 0 ~ "Buffers Positive",
    beta1 > 0 & beta2 > 0 & int_est > 0 ~ "Amplifies Positive",
    
    # Mixed main effects
    beta1 < 0 & beta2 > 0 & int_est < 0 ~ "Shifts Negative",
    beta1 > 0 & beta2 < 0 & int_est < 0 ~ "Shifts Negative",
    
    beta1 < 0 & beta2 > 0 & int_est > 0 ~ "Shifts Positive",
    beta1 > 0 & beta2 < 0 & int_est > 0 ~ "Shifts Positive",
    
    # Additive model direction
    (mode == "Additive") & beta1 < 0 & beta2 < 0 ~ "Negative",
    (mode == "Additive") & beta1 > 0 & beta2 > 0 ~ "Positive",
    (mode == "Additive") ~ "Mixed",
    
    TRUE ~ "Neutral"
  )
  
  effect_direction <- case_when(
    # --- Non-additive (interaction present) ---
    mode == "Synergistic" & beta1 < 0 & beta2 < 0 ~ "Amplifies Negative",
    mode == "Antagonistic" & beta1 < 0 & beta2 < 0 ~ "Buffers Negative",
    
    mode == "Synergistic" & beta1 > 0 & beta2 > 0 ~ "Amplifies Positive",
    mode == "Antagonistic" & beta1 > 0 & beta2 > 0 ~ "Buffers Positive",
    
    # Mixed main effects with interaction
    mode %in% c("Synergistic", "Antagonistic") &
      beta1 < 0 & beta2 > 0 & int_est < 0 ~ "Shifts Negative",
    mode %in% c("Synergistic", "Antagonistic") &
      beta1 > 0 & beta2 < 0 & int_est < 0 ~ "Shifts Negative",
    
    mode %in% c("Synergistic", "Antagonistic") &
      beta1 < 0 & beta2 > 0 & int_est > 0 ~ "Shifts Positive",
    mode %in% c("Synergistic", "Antagonistic") &
      beta1 > 0 & beta2 < 0 & int_est > 0 ~ "Shifts Positive",
    
    # --- Additive (no interaction) ---
    mode == "Additive" & beta1 < 0 & beta2 < 0 ~ "Negative (additive)",
    mode == "Additive" & beta1 > 0 & beta2 > 0 ~ "Positive (additive)",
    
    
    # Mixed effects, harmful dominates
    mode == "Additive" &
      beta1 < 0 & beta2 > 0 & abs(beta1) > abs(beta2) ~ "Shift Negative (additive)",
    
    mode == "Additive" &
      beta1 > 0 & beta2 < 0 & abs(beta1) < abs(beta2) ~ "Shift Negative (additive)",
    
    # Mixed effects, beneficial dominates
    mode == "Additive" &
      beta1 > 0 & beta2 < 0 & abs(beta1) > abs(beta2) ~ "Shift Positive (additive)",
    
    mode == "Additive" &
      beta1 < 0 & beta2 > 0 & abs(beta1) < abs(beta2) ~ "Shift Positive (additive)",
    
    # If mixed and neither dominates clearly
    mode == "Additive" ~ "Mixed (additive)",
    
    # Fallback
    TRUE ~ "Neutral"
  )
  
  # effect_direction <- case_when(
  #   beta1 < 0 & beta2 < 0 ~ "Negative",
  #   beta1 > 0 & beta2 > 0 ~ "Positive",
  #   
  #   TRUE ~ "Mixed"
  # )
  
  # ----- Output -----
  tibble(
    var1 = x,
    var2 = y,
    beta1 = beta1,
    beta2 = beta2,
    int_est = int_est,
    p_int = p_int,
    AIC_add = AIC_add,
    AIC_int = AIC_int,
    AIC_supports_add = AIC_supports_add,
    LRT_p = LRT_p,
    mode = mode,
    interaction_direction = interaction_direction,
    effect_direction = effect_direction
  )
}

