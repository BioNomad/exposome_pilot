x <- fev1_fvc_expom



x |> 
  extract_results("codebook") |> 
  filter(variable %in% num_exp_vars) |> 
  pull(category) |> 
  unique() |> 
  map(~{
    vars <- x |> 
      extract_results("codebook") |> 
      filter(category == .x) |> 
      pull(variable)
    df <- x |> 
      run_exposome_score(
        score_type = "median",
        exposure_cols = vars,
        score_column_name = paste(.x,"Score")
      ) |> 
      pivot_s()
    
    df[[paste(.x,"Score")]]
  }) |> 
  view()

x <- x |> 
  run_exposome_score(
    score_type = "median",
    exposure_cols = 
  )

# Exposure - omics --------------------
exp_omics_assoc <- fev1_fvc_expom |> 
  run_exposure_omics_association(
    exposures = num_exp_vars,
    covariates = c("gli_age","gli_sex","fis"),
    exp_name = c("CD4 T cell RNA",
                 "CD16 Monocyte RNA",
                 "Proteomics",
                 "Adduct Burden"))
### Plotting ------------------
plot_exp_class <- function(x_name){
  
  if(x_name == "CD4 T cell RNA"){
    title_lab =bquote(bold("CD4")^"+" ~ bold("T cell RNA"))
  } else if(x_name =="CD16 Monocyte RNA"){
    title_lab = bquote(bold("CD16")^"+" ~ bold("Monocyte RNA"))
  }else{
    title_lab = x_name
  }
  exp_omics_assoc |>
    extract_results("association") |> 
    pluck("exposure_omics",
          "results_df") |> 
    filter(p.value<0.05) |> 
    filter(category != "IgE Total") |> 
    dplyr::count(exposure,exp_name) |> 
    inner_join(
      exp_omics_assoc |> 
        extract_results("codebook"),
      by=c("exposure"="variable")
    ) |> 
    filter(exp_name == x_name) |> 
    group_by(category) |> 
    reframe(med_n=median(n)) |> 
    arrange(desc(med_n)) |> 
    slice_head(n=5) |> 
    ggplot(aes(
      x=med_n,
      y=reorder(category,med_n),
      fill=category
    ))+
    geom_col()+
    theme_custom()+
    scale_fill_manual(values = cat_colors)+
    theme(legend.position = "none")+
    labs(
      x="Median No. of Associations",
      y="",
      title=title_lab
    )
}
p1 <- plot_exp_class("CD4 T cell RNA")
  

p2 <- plot_exp_class("CD16 Monocyte RNA")


p3 <- plot_exp_class("Proteomics")


p4 <- plot_exp_class("Adduct Burden")
### Combined plot ---------------
library(patchwork)
((p1|p2)/(p3|p4))+plot_layout(guides = "collect")

### Median Category Rank -------------------

exp_omics_assoc |>
  extract_results("association") |> 
  pluck("exposure_omics",
        "results_df") |> 
  filter(p.value<0.05) |> 
  filter(category != "IgE Total") |> 
  dplyr::count(exposure,exp_name) |> 
  inner_join(
    exp_omics_assoc |> 
      extract_results("codebook"),
    by=c("exposure"="variable")
  ) |> 
  #filter(exp_name == "Adduct Burden") |> 
  group_by(exp_name,category) |> 
  reframe(med_n=median(n)) |> 
  arrange(desc(med_n)) |> 
  group_by(exp_name) |> 
  mutate(rank = percent_rank(med_n)) |> 
  group_by(category) |> 
  reframe(med_rank=median(rank))  |> 
  ggplot(aes(
    x=med_rank,
    y=reorder(category,med_rank),
    fill=category
  ))+
  geom_col()+
  theme_custom()+
  scale_fill_manual(values = cat_colors)+
  theme(legend.position = "none")+
  labs(
    x="Median Rank",
    y=""
  )

## Tnf Superfamily ----------------------
exp_omics_assoc |>
  extract_results("association") |> 
  pluck("exposure_omics",
        "results_df") |> 
  #filter(p.value<0.05) |> 
  filter(category != "IgE Total") |> 
  filter(exp_name == "CD4 T cell RNA") |> 
  filter(feature_id %in% c(
    msig |> 
      filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |>
      pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories |> 
      dplyr::rename(cytokine_cat=category),
    by = c("feature_id" = "gene_symbol")
  ) |> 
  filter(!cytokine_cat %in% c(
    "Misc / Less Characterized",
    "Other Chemokines"
  )) |> 
  group_by(cytokine_cat,category) |> 
  reframe(
    median_logFC = median(estimate),
    n            = dplyr::n()
  ) |> 
  filter(cytokine_cat == "TNF Superfamily") |> 
  # mutate(
  #   cytokine_cat = fct_reorder(cytokine_cat, median_logFC),
  #   direction = ifelse(median_logFC > 0, "up", "down")
  # ) |> 
  filter(n>2) |> 
  ggplot(aes(
    x    = median_logFC,
    y    = reorder(category,median_logFC),
    fill = median_logFC
  )) +
  geom_col(width = .85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey90",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Median Estimate")
  ) +
  labs(
    x     = expression("Median Estimate"),
    y     = NULL,
    subtitle = "TNF Superfamily",
    title = bquote(bold("CD4")^"+" ~ bold("T cell"))
  ) +
  theme_custom() +
  #xlim(c(-3,5.5))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    axis.text.x=element_text(angle=0,hjust=.6)
  )


exp_omics_assoc |>
  extract_results("association") |> 
  pluck("exposure_omics",
        "results_df") |> 
  #filter(p.value<0.05) |> 
  filter(category != "IgE Total") |> 
  filter(exp_name == "CD16 Monocyte RNA") |> 
  filter(feature_id %in% c(
    msig |> 
      filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |>
      pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories |> 
      dplyr::rename(cytokine_cat=category),
    by = c("feature_id" = "gene_symbol")
  ) |> 
  filter(!cytokine_cat %in% c(
    "Misc / Less Characterized",
    "Other Chemokines"
  )) |> 
  group_by(cytokine_cat,category) |> 
  reframe(
    median_logFC = median(estimate),
    n            = dplyr::n()
  ) |> 
  filter(cytokine_cat == "TNF Superfamily") |> 
  # mutate(
  #   cytokine_cat = fct_reorder(cytokine_cat, median_logFC),
  #   direction = ifelse(median_logFC > 0, "up", "down")
  # ) |> 
  filter(n>2) |> 
  ggplot(aes(
    x    = median_logFC,
    y    = reorder(category,median_logFC),
    fill = median_logFC
  )) +
  geom_col(width = .85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey90",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Median Estimate")
  ) +
  labs(
    x     = expression("Median Estimate"),
    y     = NULL,
    subtitle = "TNF Superfamily",
    title = bquote(bold("CD16")^"+" ~ bold("Monocyte"))
  ) +
  theme_custom() +
  #xlim(c(-3,5.5))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    axis.text.x=element_text(angle=0,hjust=.6)
  )

## Tgfb Superfamily ----------------------
exp_omics_assoc |>
  extract_results("association") |> 
  pluck("exposure_omics",
        "results_df") |> 
  #filter(p.value<0.05) |> 
  filter(category != "IgE Total") |> 
  filter(exp_name == "CD4 T cell RNA") |> 
  filter(feature_id %in% c(
    msig |> 
      filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |>
      pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories |> 
      dplyr::rename(cytokine_cat=category),
    by = c("feature_id" = "gene_symbol")
  ) |> 
  filter(!cytokine_cat %in% c(
    "Misc / Less Characterized",
    "Other Chemokines"
  )) |> 
  group_by(cytokine_cat,category) |> 
  reframe(
    median_logFC = median(estimate),
    n            = dplyr::n()
  ) |> 
  filter(cytokine_cat == "TGF-b / BMP / GDF") |> 
  # mutate(
  #   cytokine_cat = fct_reorder(cytokine_cat, median_logFC),
  #   direction = ifelse(median_logFC > 0, "up", "down")
  # ) |> 
  filter(n>2) |> 
  ggplot(aes(
    x    = median_logFC,
    y    = reorder(category,median_logFC),
    fill = median_logFC
  )) +
  geom_col(width = .85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey90",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Median Estimate")
  ) +
  labs(
    x     = expression("Median Estimate"),
    y     = NULL,
    subtitle = "TGF-b / BMP / GDF",
    title = bquote(bold("CD4")^"+" ~ bold("T cell"))
  ) +
  theme_custom() +
  #xlim(c(-3,5.5))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    axis.text.x=element_text(angle=0,hjust=.6)
  )


exp_omics_assoc |>
  extract_results("association") |> 
  pluck("exposure_omics",
        "results_df") |> 
  #filter(p.value<0.05) |> 
  filter(category != "IgE Total") |> 
  filter(exp_name == "CD16 Monocyte RNA") |> 
  filter(feature_id %in% c(
    msig |> 
      filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |>
      pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories |> 
      dplyr::rename(cytokine_cat=category),
    by = c("feature_id" = "gene_symbol")
  ) |> 
  filter(!cytokine_cat %in% c(
    "Misc / Less Characterized",
    "Other Chemokines"
  )) |> 
  group_by(cytokine_cat,category) |> 
  reframe(
    median_logFC = median(estimate),
    n            = dplyr::n()
  ) |> 
  filter(cytokine_cat == "TGF-b / BMP / GDF") |> 
  # mutate(
  #   cytokine_cat = fct_reorder(cytokine_cat, median_logFC),
  #   direction = ifelse(median_logFC > 0, "up", "down")
  # ) |> 
  filter(n>2) |> 
  ggplot(aes(
    x    = median_logFC,
    y    = reorder(category,median_logFC),
    fill = median_logFC
  )) +
  geom_col(width = .85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey90",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Median Estimate")
  ) +
  labs(
    x     = expression("Median Estimate"),
    y     = NULL,
    subtitle = "TGF-b / BMP / GDF",
    title = bquote(bold("CD16")^"+" ~ bold("Monocyte"))
  ) +
  theme_custom() +
  #xlim(c(-3,5.5))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    axis.text.x=element_text(angle=0,hjust=.6)
  )

## CLR + Exposure Association Testing ---------------

ps <- tidybulk::pivot_sample(decon)

cell_cols <- colnames(ps)[grepl("cell|T.reg", colnames(ps))]

cell_mat <- ps |>
  dplyr::select(all_of(cell_cols)) |>
  as.matrix()

# Zero handling (required for CLR)
cell_mat[cell_mat == 0] <- 1e-6

# CLR transform
clr_mat <- compositions::clr(cell_mat)
colnames(clr_mat) <- paste0("clr_", cell_cols)

meta <- ps |>
  dplyr::select(
    all_of(num_exp_vars),
    gli_age,
    gli_sex,
    gli_height,
    fis
  )

model_df <- bind_cols(meta, as.data.frame(clr_mat))

# Fit one model per CLR component × exposure
decon_exp_res <- map_dfr(
  colnames(clr_mat),
  function(clr_col) {
    map_dfr(
      num_exp_vars,
      function(exp_var) {
        lm(
          as.formula(paste(
            clr_col, "~",
            exp_var,
            "+ gli_age + gli_sex + fis"
          )),
          data = model_df
        ) |>
          tidy() |>
          filter(term == exp_var) |>
          mutate(
            cell = sub("^clr_", "", clr_col),
            exposure = exp_var
          )
      }
    )
  }
)

# Clean cell labels and add FDR correction
decon_exp_res <- decon_exp_res |>
  mutate(
    cell_clean = case_when(
      cell == "Memory.CD4.T.cell.TFH"             ~ "Tfh (memory CD4)",
      cell == "Memory.CD4.T.cell.Th1"             ~ "Th1 (memory CD4)",
      cell == "Memory.CD4.T.cell.Th2"             ~ "Th2 (memory CD4)",
      cell == "Memory.CD4.T.cell.Th17"            ~ "Th17 (memory CD4)",
      cell == "Memory.CD4.T.cell.Th1.Th17"        ~ "Th1/Th17 (memory CD4)",
      cell == "Terminal.effector.memory.CD4.T.cell" ~ "Terminal effector memory CD4",
      cell == "naive.CD4.T.cell"                  ~ "Naive CD4",
      cell == "T.reg"                             ~ "Treg",
      TRUE                                        ~ cell
    ),
    q.value = p.adjust(p.value, method = "BH")
  ) |> 
  left_join(
    exp_omics_assoc |> 
      extract_results("codebook"),
    by=c("exposure"="variable")
  )
#### Plot Deconvolution results -------------
##### Treg  --------------------------
treg <- decon_exp_res |> 
  group_by(category,cell_clean) |> 
  reframe(
    med_est=median(estimate),
    med_p=median(p.value)
  ) |> 
  filter(cell_clean=="Treg") |> 
  filter(category != "IgE Total") |> 
  ggplot(aes(
    x=med_est,
    y=reorder(category,med_est),
    fill=med_est
  ))+
  geom_col()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0,
    limits=c(-1,1),
    oob=scales::squish
  )+
  theme_custom()+
  xlim(c(-1.6,1.6))+
  labs(
    x="Median Estimate",
    y="",
    fill="Median Estimate",
    title = "Treg"
  )+
  theme(axis.text.x=element_text(angle=0,hjust=.6))

##### Th2  --------------------------
th2 <- decon_exp_res |> 
  group_by(category,cell_clean) |> 
  reframe(
    med_est=median(estimate),
    med_p=median(p.value)
  ) |> 
  filter(cell_clean=="Th2 (memory CD4)") |> 
  filter(category != "IgE Total") |> 
  ggplot(aes(
    x=med_est,
    y=reorder(category,med_est),
    fill=med_est
  ))+
  geom_col()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0,
    limits=c(-1,1),
    oob=scales::squish
  )+
  theme_custom()+
  xlim(c(-1.6,1.6))+
  labs(
    x="Median Estimate",
    y="",
    fill="Median Estimate",
    title = "Th2 (memory CD4)"
  )+
  theme(axis.text.x=element_text(angle=0,hjust=.6))


##### Terminal Effector Memory  --------------------------
tem <- decon_exp_res |> 
  group_by(category,cell_clean) |> 
  reframe(
    med_est=median(estimate),
    med_p=median(p.value)
  ) |> 
  filter(cell_clean=="Terminal effector memory CD4") |> 
  filter(category != "IgE Total") |> 
  ggplot(aes(
    x=med_est,
    y=reorder(category,med_est),
    fill=med_est
  ))+
  geom_col()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0,
    limits=c(-1,1),
    oob=scales::squish
  )+
  theme_custom()+
  xlim(c(-1.6,1.6))+
  labs(
    x="Median Estimate",
    y="",
    fill="Median Estimate",
    title = "Terminal effector memory CD4"
  )+
  theme(axis.text.x=element_text(angle=0,hjust=.6))

###### Combined plot ----------------

(treg|tem|th2)+plot_layout(guides="collect")
#### Tile Plot -------------------
decon_exp_res |> 
  group_by(category,cell_clean) |> 
  reframe(
    med_est=median(estimate),
    med_p=median(p.value)
  ) |> 
  ggplot(aes(
    x=reorder(category,med_est),
    y=reorder(cell_clean,med_est),
    fill=med_est
  ))+
  geom_tile()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0,
    limits=c(-.5,.5),
    oob=scales::squish
  )+
  theme_custom()



  ggplot(aes(x = estimate, 
             y = fct_reorder(cell_clean,estimate),
             color = estimate,
             fill = estimate)) +
  geom_point(shape = 18, 
             size = 4, 
             stroke = 0.8) +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error),
                 height = 0.3) +
  scale_color_gradient2(
    low="#5399b0",
    mid="grey85",
    high="#9B6981FF",
    midpoint = 0
  )+
  scale_fill_gradient2(
    low="#5399b0",
    mid="grey85",
    high="#9B6981FF",
    midpoint = 0
  )+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = NULL,
    x = paste("Effect size",
              "(CLR log-ratio, 95% CI)",
              sep = "\n"),
    y = "",
    fill ="Effect size",
    color ="Effect size"
  ) +
  theme_bw(base_size = 12) +
  theme_custom()+
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle=0,
                                    face="bold.italic"),
        axis.text.y = element_text(face = "plain"),
        legend.position = "none")