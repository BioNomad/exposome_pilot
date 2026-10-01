## ---------------
combo_assoc_df <- fev1_fvc_expom |> 
  extract_results(result = "association") |> 
  pluck("assoc_exposures",
        "results_df") |> 
  dplyr::select(term,
                estimate,
                std.error,
                p.value,
                category) |> 
  left_join(sens_df |> 
              dplyr::select(exposure,
                            median_effect,
                            sign_consistency),
            by=c("term"="exposure")) |> 
  left_join(broom::tidy(enet_res$glmnet.fit) |>
              filter(lambda == enet_res$lambda.min,
                     term != "(Intercept)") |> 
              dplyr::select(term,
                            enet_estimate=estimate)) |> 
  mutate(
    sig=p.value<0.05,
    term = case_when(
      sig == TRUE ~ paste0("* ",term),
      .default =term),
    term = forcats::fct_reorder(term, estimate),
    category = forcats::fct_reorder(
      category, estimate, .fun = median, .desc = TRUE)
  )

p_assoc <-
  combo_assoc_df |>
  ggplot(aes(x = estimate,
             y = term,
             color = estimate,
             fill = estimate)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.8) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error),
    height = 0.15
  ) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  scale_color_gradient2(
    low = "blue4", mid = "grey85", high = "red4", midpoint = 0
  ) +
  scale_fill_gradient2(
    low = "blue4", mid = "grey85", high = "red4", midpoint = 0
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size", y = NULL,color="Effect Size") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold.italic"),
        legend.position = "bottom")

p_stability <-
  combo_assoc_df |>
  ggplot(aes(x = sign_consistency,
             y = term,
             color = median_effect)) +
  geom_segment(
    aes(x = 0, xend = sign_consistency,
        yend = term),
    linewidth = 0.8
  ) +
  geom_point(size = 2.5) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  scale_color_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  labs(x = "Effect consistency", y = NULL, color = "Median Effect") +
  theme_bw(base_size = 13)+
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold.italic"),
        legend.position = "bottom")

p_enet <-
  combo_assoc_df |>
  ggplot(aes(x = enet_estimate,
             y = term,
             color = enet_estimate)) +
  geom_segment(
    aes(x = 0, xend = enet_estimate,
        yend = term),
    linewidth = 0.8
  ) +
  geom_point(size = 2.5) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  labs(x = "Elastic-net coefficient", y = NULL,,color="Effect Size") +
  theme_bw(base_size = 13)+
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold.italic"),
        legend.position = "bottom")

p_assoc | p_stability | p_enet


## --------------------------------------------------------------------------------------------
combo_assoc_df |>
  ggplot(aes(x = enet_estimate,
             y = term,
             color = enet_estimate)) +
  geom_segment(
    aes(x = 0, xend = enet_estimate,
        yend = term),
    linewidth = 0.8
  ) +
  geom_point(size = 2.5) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  labs(x = "Elastic-net coefficient", y = NULL,,color="Effect Size") +
  theme_bw(base_size = 13)+
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold.italic"),
        legend.position = "right")


## -------------
combo_assoc_df <- fev1_fvc_expom |> 
  extract_results(result = "association") |> 
  pluck("assoc_exposures",
        "results_df") |> 
  dplyr::select(term,
                estimate,
                std.error,
                p.value,
                category) |> 
  left_join(sens_df |> 
              filter(sign_consistency>.5) |> 
              dplyr::select(exposure,
                            median_effect,
                            sign_consistency),
            by=c("term"="exposure")) |> 
  left_join(broom::tidy(enet_res$glmnet.fit) |>
              filter(lambda == enet_res$lambda.min,
                     term != "(Intercept)") |> 
              dplyr::select(term,
                            enet_estimate=estimate)) |> 
  mutate(
    sig=p.value<0.05,
    term = case_when(
      sig == TRUE ~ paste0("* ",term),
      .default =term),
    term = forcats::fct_reorder(term, estimate),
    category = forcats::fct_reorder(
      category, estimate, .fun = median, .desc = TRUE)
  )

p_assoc <-
  combo_assoc_df |>
  ggplot(aes(x = estimate,
             y = term,
             color = estimate,
             fill = estimate)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.8) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error),
    height = 0.15
  ) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  scale_color_gradient2(
    low = "blue4", mid = "grey85", high = "red4", midpoint = 0
  ) +
  scale_fill_gradient2(
    low = "blue4", mid = "grey85", high = "red4", midpoint = 0
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size", y = NULL,color="Effect Size",fill="Effect Size") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background.y = element_blank())

p_stability <-
  combo_assoc_df |>
  ggplot(aes(x = sign_consistency,
             y = term,
             color = median_effect)) +
  geom_segment(
    aes(x = 0, xend = sign_consistency,
        yend = term),
    linewidth = 0.8
  ) +
  geom_point(size = 2.5) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  scale_color_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish
  ) +scale_x_continuous(breaks = c(0,0.5,1))+
  labs(x = "Effect consistency", y = NULL, color = "Median Effect") +
  theme_bw(base_size = 13)+
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),)

p_enet <-
  combo_assoc_df |>
  ggplot(aes(x = enet_estimate,
             y = term,
             color = enet_estimate)) +
  geom_segment(
    aes(x = 0, xend = enet_estimate,
        yend = term),
    linewidth = 0.8
  ) +
  geom_point(size = 2.5) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,
    limits = c(-0.05, 0.05),
    oob = scales::squish
  ) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "grey85",
    high = "red4",
    midpoint = 0,na.value = "white",
    limits = c(-0.05, 0.05),
    oob = scales::squish
  ) +
  labs(x = "Elastic-net coefficient", y = NULL,fill="Effect Size",color="Effect Size") +
  theme_bw(base_size = 13)+
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold.italic"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

(p_assoc | p_stability | p_enet) +
  patchwork::plot_layout(widths = c(3, 3, 3))


## ------------------
combo_assoc_df <- fev1_fvc_expom |> 
  extract_results(result = "association") |> 
  pluck("assoc_exposures",
        "results_df") |> 
  dplyr::select(term,
                estimate,
                std.error,
                p.value,
                category) |> 
  left_join(sens_df |> 
              #filter(sign_consistency>.5) |> 
              dplyr::select(exposure,
                            median_effect,
                            sign_consistency),
            by=c("term"="exposure")) |> 
  left_join(broom::tidy(enet_res$glmnet.fit) |>
              filter(lambda == enet_res$lambda.min,
                     term != "(Intercept)") |> 
              dplyr::select(term,
                            enet_estimate=estimate)) |> 
  mutate(
    sig=p.value<0.05,
    term = case_when(
      sig == TRUE ~ paste0("* ",term),
      .default =term),
    term = forcats::fct_reorder(term, estimate),
    category = forcats::fct_reorder(
      category, estimate, .fun = median, .desc = TRUE)
  )

p_assoc <-
  combo_assoc_df |>
  ggplot(aes(x = estimate,
             y = term,
             color = estimate,
             fill = estimate)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.8) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error),
    height = 0.15
  ) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  scale_color_gradient2(
    low = "blue4", mid = "grey85", high = "red4", midpoint = 0
  ) +
  scale_fill_gradient2(
    low = "blue4", mid = "grey85", high = "red4", midpoint = 0
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size", y = NULL,color="Effect Size",fill="Effect Size") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        strip.text.y = element_blank(),
        strip.background.y = element_blank())

p_stability <-
  combo_assoc_df |>
  ggplot(aes(x = 1,
             y = term,
             color = sign_consistency)) +
  geom_point(aes(size=sign_consistency)) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  scale_color_gradientn(
    colors = c(
      "#3a1c71",
      "#d76d77",
      "#ffaf7b"
    ),
    na.value = "grey80")+
  scale_x_continuous(breaks = c(0,0.5,1))+
  labs(x = "Effect consistency", y = NULL, color = "Effect Consistency") +
  theme_void(base_size = 13)+guides(size="none")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),)

p_enet <-
  combo_assoc_df |>
  ggplot(aes(x = 1,
             y = term,
             color = enet_estimate)) +
  geom_point(aes(size = enet_estimate)) +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space = "free",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(assoc_df$category)]
      )
    )
  ) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(
    low = "#5D26C1",
    mid = "grey85",
    high = "#4BC0C8",
    midpoint = 0,
    limits = c(-0.05, 0.05),
    oob = scales::squish
  ) +
  labs(x = "Elastic-net coefficient", y = NULL,color="Enet Effect Size") +
  theme_void(base_size = 13)+guides(size="none")+
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold.italic"),
        legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

(p_assoc | p_stability | p_enet) +
  patchwork::plot_layout(widths = c(3, 1, 1),guides = "collect")