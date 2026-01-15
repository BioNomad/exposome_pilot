summarise_set <- function(df, genes, cell_pattern, set_label) {
  df |>
    filter(
      feature_map %in% genes,
      grepl(cell_pattern, exp_name)
    ) |>
    group_by(exposure) |>
    reframe(
      med_signed = median(signed_logFC, na.rm = TRUE),
      frac_expected_dir = mean(signed_logFC > 0, na.rm = TRUE),
      n_genes = n_distinct(feature_map)
    ) |>
    filter(
      n_genes > 1,
      frac_expected_dir >= 0.5,
      med_signed > 0
    ) |>
    mutate(set = set_label)
}

senescence_gate <- function(df, cell_pattern) {
  df |>
    filter(
      feature_map %in% c("CDKN2A", "CDKN1A"),
      grepl(cell_pattern, exp_name)
    ) |>
    group_by(exposure) |>
    reframe(gate_signal = median(logFC, na.rm = TRUE)) |>
    filter(gate_signal > 0)
}
