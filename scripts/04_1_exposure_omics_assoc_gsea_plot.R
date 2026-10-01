
p1 <- exp_gsea_res_sig_clean |> 
  filter(exp_name == "CD4 T cell RNA") |> 
  group_by(term_name) |>
  reframe(
    n_dir      = n_distinct(sign(NES)),
    n_exposure = n_distinct(exposure),
    med_nes    = median(NES),
    exp_name   = exp_name
  ) |>
  filter(n_dir==1) |> 
  arrange(med_nes)  |> 
  distinct() |> 
  slice_max(abs(med_nes), n = 1, by = term_name) |>
  mutate(rank = row_number()) |> 
  (\(df) {
    label_df <- df |>
      filter(exp_name == "CD4 T cell RNA", term_name %in%  c(
        "endoplasmic reticulum unfolded protein response",
        "cellular response to hypoxia",
        "regulation of mRNA splicing, via spliceosome",
        "cell chemotaxis",
        "mitochondrial inner membrane",
        "cytoplasmic stress granule",
        "canonical Wnt signaling pathway",
        "transcription coregulator activity",
        "cell-cell adhesion",
        "alpha-beta T cell receptor complex",
        "T cell activation"
      )) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(3, 12, length.out = n()))
    
    df |>
      mutate(
        fill_nes = if_else(exp_name == "CD4 T cell RNA", med_nes, NA_real_),
        facet_y  = "CD4 T cell RNA"
      ) |>
      ggplot(aes(x = rank, y = facet_y, fill = fill_nes)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos, y = 1.5, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos + 0.5, y = y_label, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_label(
        data          = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes   = FALSE,
        hjust         = 0,
        size          = 4.5,
        color         = "grey20",
        fill          = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits   = c(-2, 2),
        oob      = scales::squish,
        na.value = "grey92"
      ) +
      coord_cartesian(clip = "off") +
      theme_custom() +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "CD4 T cell Isoforms"       = expression(CD4^"+" ~ "T cell Isoforms"),
        "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
        "Proteomics" ~ "Proteomics",
        "Adductomics" ~ "Adductomics",
        "Adduct Burden" ~ "Adduct Burden"
      ))+
      theme(
        axis.text.x           = element_blank(),
        axis.line.y           = element_blank(),
        axis.line.x           = element_blank(),
        axis.ticks.y          = element_blank(),
        axis.ticks.x          = element_blank(),
        legend.position       = "bottom",
        legend.title.position = "top",
        plot.margin           = margin(40, 300, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "Median NES")
  })()

p2 <- exp_gsea_res_sig_clean |> 
  filter(exp_name == "CD16 Monocyte RNA") |> 
  group_by(term_name) |>
  reframe(
    n_dir      = n_distinct(sign(NES)),
    n_exposure = n_distinct(exposure),
    med_nes    = median(NES),
    exp_name   = exp_name
  ) |>
  filter(n_dir==1) |> 
  arrange(med_nes)  |> 
  distinct() |> 
  slice_max(abs(med_nes), n = 1, by = term_name) |>
  mutate(rank = row_number()) |>
  (\(df) {
    label_df <- df |>
      filter(exp_name == "CD16 Monocyte RNA", term_name %in% c(
        "defense response to virus",
        "signaling receptor binding",
        "toll-like receptor 9 signaling pathway",
        "positive regulation of type II interferon production"
      )) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(2, 5, length.out = n()))
    
    df |>
      mutate(
        fill_nes = if_else(exp_name == "CD16 Monocyte RNA", med_nes, NA_real_),
        facet_y  = "CD16 Monocyte RNA"
      ) |>
      ggplot(aes(x = rank, y = facet_y, fill = fill_nes)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos, y = 1.5, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos + 0.5, y = y_label, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_label(
        data          = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes   = FALSE,
        hjust         = 0,
        size          = 4.5,
        color         = "grey20",
        fill          = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits   = c(-2, 2),
        oob      = scales::squish,
        na.value = "grey92"
      ) +
      coord_cartesian(clip = "off") +
      theme_custom() +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "CD4 T cell Isoforms"       = expression(CD4^"+" ~ "T cell Isoforms"),
        "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
        "Proteomics" ~ "Proteomics",
        "Adductomics" ~ "Adductomics",
        "Adduct Burden" ~ "Adduct Burden"
      ))+
      theme(
        axis.text.x           = element_blank(),
        axis.line.y           = element_blank(),
        axis.line.x           = element_blank(),
        axis.ticks.y          = element_blank(),
        axis.ticks.x          = element_blank(),
        legend.position       = "bottom",
        legend.title.position = "top",
        plot.margin           = margin(40, 200, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "Median NES")
  })()

p3 <- exp_gsea_res_sig_clean |> 
  filter(exp_name == "Proteomics") |> 
  group_by(term_name) |>
  reframe(
    n_dir      = n_distinct(sign(NES)),
    n_exposure = n_distinct(exposure),
    med_nes    = median(NES),
    exp_name   = exp_name
  ) |>
  filter(n_dir==1) |> 
  arrange(med_nes)  |> 
  distinct() |> 
  slice_max(abs(med_nes), n = 1, by = term_name) |>
  mutate(rank = row_number()) |>
  (\(df) {
    label_df <- df |>
      filter(exp_name == "Proteomics", term_name %in% c(
        "secretory granule lumen",
        "cell-cell adhesion",
        "complement activation, classical pathway",
        "immunoglobulin receptor binding"
      )) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(2, 5, length.out = n()))
    
    df |>
      mutate(
        fill_nes = if_else(exp_name == "Proteomics", med_nes, NA_real_),
        facet_y  = "Proteomics"
      ) |>
      ggplot(aes(x = rank, y = facet_y, fill = fill_nes)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos, y = 1.5, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos + 0.5, y = y_label, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_label(
        data          = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes   = FALSE,
        hjust         = 0,
        size          = 4.5,
        color         = "grey20",
        fill          = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits   = c(-2, 2),
        oob      = scales::squish,
        na.value = "grey92"
      ) +
      coord_cartesian(clip = "off") +
      theme_custom() +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "CD4 T cell Isoforms"       = expression(CD4^"+" ~ "T cell Isoforms"),
        "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
        "Proteomics" ~ "Proteomics",
        "Adductomics" ~ "Adductomics",
        "Adduct Burden" ~ "Adduct Burden"
      ))+
      theme(
        axis.text.x           = element_blank(),
        axis.line.y           = element_blank(),
        axis.line.x           = element_blank(),
        axis.ticks.y          = element_blank(),
        axis.ticks.x          = element_blank(),
        legend.position       = "bottom",
        legend.title.position = "top",
        plot.margin           = margin(40, 200, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "Median NES")
  })()

p4 <- exp_gsea_res_sig_clean |> 
  filter(exp_name == "Adduct Burden") |> 
  group_by(term_name) |>
  reframe(
    n_dir      = n_distinct(sign(NES)),
    n_exposure = n_distinct(exposure),
    med_nes    = median(NES),
    exp_name   = exp_name
  ) |>
  filter(n_dir==1) |> 
  distinct() |> 
  arrange(med_nes)  |> 
  slice_max(abs(med_nes), n = 1, by = term_name) |>
  mutate(rank = row_number()) |>
  (\(df) {
    label_df <- df |>
      filter(exp_name == "Adduct Burden", 
             term_name %in% c(
        "extracellular region",
        "blood microparticle",
        "RNA binding"
      )) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(2, 5, length.out = n()))
    
    df |>
      mutate(
        fill_nes = if_else(exp_name == "Adduct Burden", med_nes, NA_real_),
        facet_y  = "Adduct Burden"
      ) |>
      ggplot(aes(x = rank, y = facet_y, fill = fill_nes)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos, y = 1.5, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_segment(
        data        = label_df,
        aes(x = x_pos, xend = x_pos + 0.5, y = y_label, yend = y_label),
        inherit.aes = FALSE,
        color       = "grey40",
        linewidth   = 0.3
      ) +
      geom_label(
        data          = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes   = FALSE,
        hjust         = 0,
        size          = 4.5,
        color         = "grey20",
        fill          = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits   = c(-2, 2),
        oob      = scales::squish,
        na.value = "grey92"
      ) +
      coord_cartesian(clip = "off") +
      theme_custom() +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "CD4 T cell Isoforms"       = expression(CD4^"+" ~ "T cell Isoforms"),
        "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
        "Proteomics" ~ "Proteomics",
        "Adductomics" ~ "Adductomics",
        "Adduct Burden" ~ "Adduct Burden"
      ))+
      theme(
        axis.text.x           = element_blank(),
        axis.line.y           = element_blank(),
        axis.line.x           = element_blank(),
        axis.ticks.y          = element_blank(),
        axis.ticks.x          = element_blank(),
        legend.position       = "bottom",
        legend.title.position = "top",
        plot.margin           = margin(40, 200, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "Median NES")
  })()

library(patchwork)
(p1 / p2 / p3 / p4) +
  plot_layout(guides = "collect", heights = c(3, 1, 1, 1))&
  theme(legend.position = "bottom")

# old col plot --------
# make_bp <- function(exp) {
#   consistent_terms <- exp_gsea_res_sig_clean |>
#     filter(exp_name == exp) |>
#     group_by(term_name) |>
#     reframe(
#       n_dir    = n_distinct(sign(NES)),
#       exp_name = exp_name
#     ) |>
#     filter(n_dir == 1) |>
#     distinct()
#   
#   plot_data <- exp_gsea_res_sig_clean |>
#     filter(exp_name == exp) |>
#     filter(term_name %in% consistent_terms$term_name) |>
#     filter(exposure %in% num_exp_vars)
#   
#   n_exp <- plot_data |>
#     dplyr::count(exposure) |>
#     nrow()
#   
#   plot_data |>
#     group_by(category, exposure) |>
#     summarise(n_pathways = n_distinct(term_name), .groups = "drop") |>
#     group_by(category) |>
#     summarise(med_pathways = median(n_pathways), .groups = "drop") |>
#     right_join(
#       aw_cb |>
#         filter(variable %in% num_exp_vars) |>
#         dplyr::select(category) |>
#         distinct(),
#       by = "category"
#     ) |>
#     replace_na(list(med_pathways = 0)) |>
#     ggplot(aes(
#       x    = category,
#       y    = med_pathways,
#       fill = category
#     )) +
#     geom_col() +
#     theme_custom() +
#     scale_fill_manual(values = cat_colors) +
#     theme(
#       axis.text.x           = element_blank(),
#       axis.ticks.x          = element_blank(),
#       axis.line.x           = element_blank(),
#       legend.position       = "bottom",
#       legend.title.position = "top"
#     ) +
#     #guides(fill = "none") +
#     ylim(c(0,115))+
#     labs(
#       x       = "",
#       y       = "Median pathways affected",
#       fill    = "Category",
#       title   = exp,
#       caption = paste("No. of exposures =", n_exp)
#     )
# }
# 
# bp1 <- make_bp("CD4 T cell RNA")
# bp2 <- make_bp("CD16 Monocyte RNA")
# bp3 <- make_bp("Proteomics")
# bp4 <- make_bp("Adduct Burden")
# 

# New ---------
make_bp <- function(exp) {
  consistent_terms <- exp_gsea_res_sig_clean |>
    filter(exp_name == exp) |>
    group_by(term_name) |>
    reframe(
      n_dir      = n_distinct(sign(NES)),
      n_exposure = n_distinct(exposure),
      med_nes    = median(NES),
      exp_name   = exp_name
    ) |>
    filter(n_dir == 1) |>
    distinct()
  
  n_exp <- exp_gsea_res_sig_clean |>
    filter(exp_name == exp) |>
    filter(term_name %in% consistent_terms$term_name) |>
    dplyr::count(exposure) |>
    nrow()
  
  med_pathways <- exp_gsea_res_sig_clean |>
    filter(exp_name == exp) |>
    filter(term_name %in% consistent_terms$term_name) |>
    filter(exposure %in% num_exp_vars) |>
    group_by(category, exposure) |>
    summarise(n_pathways = n_distinct(term_name), .groups = "drop") |>
    group_by(category) |>
    summarise(med_pathways = median(n_pathways), .groups = "drop")
  
  exp_gsea_res_sig_clean |>
    filter(exp_name == exp) |>
    filter(term_name %in% consistent_terms$term_name) |>
    distinct(category) |>
    mutate(present = "present") |>
    right_join(
      aw_cb |>
        filter(variable %in% num_exp_vars) |>
        dplyr::select(category),
      by = "category"
    ) |>
    distinct() |>
    left_join(med_pathways, by = "category") |>
    arrange(category) |>
    ggplot(aes(
      x    = category,
      y    = 1,
      fill = category
    )) +
    geom_tile(aes(alpha = present)) +
    geom_point(
      aes(size = med_pathways),
      shape = 21,fill="grey95",
      colour = "black"
    ) +
    theme_custom() +
    scale_fill_manual(values = cat_colors) +
    scale_alpha_manual(
      values   = c("present" = 1),
      na.value = 0.0000001
    ) +
    scale_size_continuous(range = c(1, 8)) +
    theme(
      axis.text.x           = element_blank(),
      axis.text.y           = element_blank(),
      axis.ticks            = element_blank(),
      axis.line             = element_blank(),
      legend.position       = "bottom",
      legend.title.position = "top"
    ) +
    guides(alpha = "none" ) +
    labs(
      x       = "",
      y       = "",
      size    = "Median pathways",
      title   = exp,
      fill = "Category",
      caption = paste("No. of exposures =", n_exp)
    )
}

bp1 <- make_bp("CD4 T cell RNA")
bp2 <- make_bp("CD16 Monocyte RNA")
bp3 <- make_bp("Proteomics")
bp4 <- make_bp("Adduct Burden")

