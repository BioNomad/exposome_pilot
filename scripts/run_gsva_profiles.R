run_gsva_profiles <- function(
  gs_names,
  expom,
  exp_name,
  outcome,
  covariates
  
  ){
  gois <- msig |> 
    filter(gs_name %in% gs_names) |> 
    (\(df) split(df$gene_symbol,df$gs_name))()
  
  profiles <- get_gsva_scores(
    tidyexposomics:::.update_assay_colData(
      normalize_rna_assay(expom,exp_name),
      exp_name
    ),
    gene_sets = gois,
    gene_symbol_col = "feature_id",
    log_trans = F
  )
  
  lm_res <- map(names(gois),~{
    lm(
      as.formula(paste0(
        outcome,
        " ~ ",
        #"pftfev1fvc_actual ~ ",
        .x,
        " + ",
        paste(covariates,collapse = " + ")
        #" + gli_age + gli_sex + fis"
      )),
      data = profiles
    ) |> 
      broom::tidy() 
  }) |> 
    bind_rows() |> 
    filter(term %in% names(gois)) |> 
    mutate(exp_name = exp_name)
  
  res_lst <- list(
    profiles=profiles,
    lm_res=lm_res
  )
}


