get_viper <- function(
    deg_df,
    lvls = c("A", "B", "C", "D"),
    meta_analysis = FALSE
) {
  library(dorothea)
  library(viper)
  
  dorothea_hs <- dorothea::dorothea_hs |> filter(confidence %in% lvls)
  viper_regulons <- df2regulon(dorothea_hs)
  
  if (meta_analysis) {
    signature       <- with(deg_df, qnorm(pooled_p / 2, lower.tail = FALSE) * sign(pooled_logFC))
    names(signature) <- deg_df$feature_id
  } else {
    signature       <- with(deg_df, qnorm(pvalue / 2, lower.tail = FALSE) * sign(logfc))
    names(signature) <- deg_df$feature_id
  }
  
  mrs <- msviper(signature, viper_regulons)
  mrs <- ledge(mrs)
  
  ledge_tbl <- mrs$ledge |>
    map(\(l) paste(l, collapse = ",")) |>
    tibble::enframe(name = "tf", value = "ledge_genes_full")
  
  summary(mrs, length(mrs$es$nes)) |>
    tibble::rownames_to_column("tf") |>
    tibble::as_tibble() |>
    janitor::clean_names() |>
    dplyr::select(-ledge, -regulon) |>
    left_join(ledge_tbl, by = "tf")
}