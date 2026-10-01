## Load Libraries -------------------------
library(tidyverse)
library(janitor)
library(ggpubr)
library(tidyexposomics)
library(MultiAssayExperiment)
library(biostats)
library(glmnet)
library(gtsummary)

source("./scripts/bin/internals.R")
source("./scripts/bin/labels_clean.R")

c(
  high = "#9B6981FF",
  mid = "white",
  low = "#5399b0"
)

## Load Data --------------------------------
# MultiAssayExperiment Object
expom <- readRDS("./results/qc/expom.rds")

# Codebook - Filtered to include variables
# retained in MultiAssayExperiment Object
aw_cb <- readRDS("./results/input_data/aw_cb.rds") |> 
  filter(variable %in% colnames(colData(expom)))

# Exposure variables
exp_vars <- readRDS("./results/qc/exp_vars.rds")
num_exp_vars <- readRDS("./results/qc/num_exp_vars.rds")

# All variables
vars <- aw_cb$variable |> 
  (\(chr)chr[chr %in% colnames(colData(expom))])()

# Setting colors for exposure categories
cat_colors <- scales::alpha(get_palette("aaas",13), 0.5) |> 
  (\(colors){
    names(colors) <- aw_cb |> 
      filter(variable %in% exp_vars) |> 
      pull(category) |> 
      unique()
    colors
  })()


## Sample Clustering -----------------


# Sample clustering
a <- expom |>
    run_cluster_samples(
        exposure_cols = num_exp_vars,
        clustering_approach = "elbow",
        action = "add"
    )

a |>
    plot_sample_clusters(
        exposure_cols = num_exp_vars
    )


cor_mat <- expom |> 
  pivot_sample() |> 
  dplyr::select(num_exp_vars, .sample) |> 
  column_to_rownames(".sample") |> 
  t() |> 
  cor()

# get clusters from a first pass
clusters <- cor_mat |> 
  pheatmap::pheatmap(silent = TRUE) |> 
  _$tree_col |> 
  cutree(k = 2) |> 
  as.data.frame() |> 
  setNames("Cluster") |> 
  mutate(Cluster = factor(paste("Cluster",Cluster)))

# final plot with annotation
cor_mat |> 
  pheatmap::pheatmap(
    color = colorRampPalette(c("#1D2671", "white", "#C33764"))(50),
    annotation_col = clusters,
    annotation_row = clusters,
    annotation_colors = list(
      Cluster = c("Cluster 1" = "#bbd2c5", 
                  "Cluster 2" = "#191654")
    )
  )

## Why are there two clusters? -------------


ph <- expom |> 
  pivot_sample() |> 
  dplyr::select(num_exp_vars, .sample) |> 
  column_to_rownames(".sample") |> 
  t() |> 
  cor() |> 
  pheatmap::pheatmap()


cluster_df <- cutree(ph$tree_col, k = 2) |> 
  as.data.frame() |> 
  setNames("cluster") |> 
  rownames_to_column(".sample")

expom |> 
  pivot_sample() |> 
  dplyr::select(num_exp_vars, .sample) |> 
  left_join(clusters |> 
              rownames_to_column(".sample"),
            by = ".sample") |> 
  pivot_longer(num_exp_vars, 
               names_to = "variable", 
               values_to = "value") |> 
  group_by(variable, Cluster) |> 
  reframe(mean = mean(value, na.rm = TRUE)) |> 
  pivot_wider(names_from = Cluster, 
              values_from = mean) |> 
  mutate(diff = `Cluster 2` - `Cluster 1`) |> 
  arrange(desc(abs(diff))) |> 
  left_join(
    extract_results(expom,result="codebook"),
    by="variable"
  ) |> 
  ggplot(aes(
    x=1,
    y=reorder(clean_name,diff),
    fill=diff
  ))+
  geom_tile()+
  theme_custom()+
  scale_fill_gradient2(
    high = "#9B6981FF",
    mid = "grey85",
    low = "#5399b0",
    midpoint = 0,
  )+
  labs(
    x="",
    y="",
    fill = expression(~Delta*" Cluster 2 - Cluster 1")
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank()
  )

## Correlated Exposures -----------------

num_vars <- expom@colData |>
  as_tibble() |>
  select_if(is.numeric) |> 
  (\(df){df[,colnames(df) %in% vars]})() |> 
  colnames()

expom <- expom |>
    run_correlation(
        feature_type = "exposures",
        action = "add",
        exposure_cols = num_vars,
        correlation_cutoff = 0.3
    )

expom |>
    plot_circos_correlation(
        feature_type = "exposures",
        corr_threshold = 0.5,
        exposure_cols = num_vars,
        annotation_colors = get_palette("jco",k=42) |> 
          (\(colors){
            names(colors) <- unique(aw_cb$category)
            colors
          })()
        
    )


## Correlation Circos ---------------------------
cors <- expom |>
    run_correlation(
        feature_type = "exposures",
        action = "get",
        exposure_cols = exp_vars,
        correlation_method = "spearman",
        correlation_cutoff = 0,
        pval_cutoff = 1
    )


cor_g <- cors |>
    filter(
        p.value < 0.05,
        abs(correlation) > 0.3
    ) |>
    tidygraph::as_tbl_graph(directed = FALSE) |>
    tidygraph::activate(nodes) |>
    mutate(
        community = tidygraph::group_louvain(
            weights = abs(correlation)
        ),
        strength = tidygraph::centrality_degree(
            weights = abs(correlation)
        ),
        top_lab =if_else(strength > 10, name, NA_character_)
    ) |> 
    left_join(aw_cb,
              by=c("name"="variable")) |> 
  tidygraph::rename(Category=category)

cor_g <- cor_g |> 
    group_by(Category) |> 
    arrange(desc(strength), .by_group = TRUE) |>
    mutate(
        top_lab = if_else(name %in% c(
            cor_g |>
                as.data.frame() |> 
                group_by(Category) |> 
                arrange(desc(strength),
                        .by_group = TRUE) |> 
                tidygraph::slice_head(n=2) |>
                pull(name)
        ), name, NA_character_)
    ) |>
    ungroup()



cor_g |>
    ggraph::ggraph(
        layout = "fr"
    ) +
    ggraph::geom_edge_link(
        aes(
            color = correlation,
            width = abs(correlation)*5
        )
    ) +
    ggraph::geom_node_point(aes(fill=Category,size=strength),
                            shape = 21,
                            #size = 4,
                            alpha = 0.8
    ) +
    scale_fill_manual(values=cat_colors)+
    ggraph::geom_node_label(
        aes(
            label = top_lab
        ),
        fontface = "bold.italic",repel = T,
        na.rm = TRUE
    ) +
    ggraph::theme_graph() +
    ggraph::scale_edge_color_gradient2(
        high = "#C33764",
        mid  = "white",
        low  = "#1D2671",
        midpoint = 0.1,
        guide = ggraph::guide_edge_colorbar()
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(size = 4)
        ),size = "none",
        edge_width = "none"
    )+
    ggplot2::labs(edge_color = "Correlation")



## Individual Exposure Correlations -----------------
cors |>
    select(var1, var2, correlation) |>
    
    # --- Build correlation matrix ---
    pivot_wider(names_from = var2, values_from = correlation) |>
    column_to_rownames("var1") |>
    (\(mat){ mat[is.na(mat)] <- 0; mat })() |>
    
    # --- Wrap into pheatmap call ---
    (\(mat){
        
        # Row annotation (aw_cb$type mapped to matching variables)
        row_ann <- aw_cb |>
            select(variable, category) |>
            filter(variable %in% rownames(mat)) |>
            distinct() |>
            column_to_rownames("variable")
        
        # Column annotation
        col_ann <- aw_cb |>
            select(variable, category) |>
            filter(variable %in% colnames(mat)) |>
            distinct() |>
            column_to_rownames("variable")
        
        # Colors for main heatmap
        heat_cols <- colorRampPalette(c("#1d2671","white","#c33764"))(100)
        
        # Category colors for annotations
        cats <- unique(c(row_ann$category, col_ann$category))
        ann_cols <- list(
            # category = setNames(
            #     scales::alpha(get_palette("aaas",length(cats)), 0.5),
            #     cats
            # )
            category = cat_colors
        )
        
        pheatmap::pheatmap(
            mat,
            color = heat_cols,
            breaks = seq(-1, 1, length.out = 101),
            border_color = NA,
            annotation_row = row_ann,
            annotation_col = col_ann,
            annotation_colors = ann_cols,
            fontsize = 9
        )
    })()



## Summarize Exposure Correlations --------------
cor_sum_df <- cors |>
  select(var1, var2, correlation) |>

  # --- Bring categories into long form first ---
  left_join(aw_cb |> select(variable, category),
            by = c("var1" = "variable")) |>
  dplyr::rename(cat1 = category) |>
  left_join(aw_cb |> select(variable, category),
            by = c("var2" = "variable")) |>
  dplyr::rename(cat2 = category) |>

  # --- Drop rows without category info ---
  filter(!is.na(cat1), !is.na(cat2)) |>

  # --- Label within vs. between correlation types ---
  mutate(
    corr_type = if_else(cat1 == cat2, "Within-category", "Between-category"),
    pair_category = if_else(cat1 == cat2, cat1, "Between")
  ) 


cor_sum_df |> 
  mutate(combo=paste(cat1,cat2)) |> 
  group_by(combo) |> 
  reframe(med_cor=median(correlation),cat1,cat2) |>
  distinct() |>
  dplyr::select(cat1,cat2,med_cor) |> 
  arrange(med_cor) |> 
    tidygraph::as_tbl_graph() |> 
    ggraph::ggraph(
        layout = "linear",
        circular = TRUE
    ) +
    ggraph::geom_edge_arc(aes(
        color = med_cor,
        width = abs(med_cor)*10
    )) +
    ggraph::geom_node_point(shape = 21,
                            size = 4,
                            alpha = 0.8
    ) +
    ggraph::geom_node_text(
        aes(
            label = name,
            x = 1.1 * x,
            y = 1.1 * y,
            angle = ggraph::node_angle(x, y)
        ),
        hjust = "outward",
        fontface = "bold.italic",
        check_overlap = TRUE
    ) +
    ggraph::theme_graph() +
    ggraph::scale_edge_color_gradient2(
        high = "#C33764",
        mid = "white",
        low = "#1D2671",
        midpoint = 0.1,
        guide = ggraph::guide_edge_colorbar()
    ) +
    ggplot2::coord_fixed(xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::guides(edge_width = "none") +
    ggplot2::labs(
        edge_color = "Median Correlation"
    )

## Exposure  Correlation Heatmap --------------
cor_sum_df |>
  mutate(combo = paste(cat1, cat2)) |>
  group_by(combo) |>
  reframe(med_cor = median(correlation), cat1, cat2) |>
  distinct() |>
  dplyr::select(cat1, cat2, med_cor) |>
  (\(df) {
    mat <- df |>
      tidyr::pivot_wider(names_from = cat2, values_from = med_cor, values_fill = 0) |>
      tibble::column_to_rownames("cat1") |>
      as.matrix()
    ord <- hclust(dist(mat))$order
    lvls <- rownames(mat)[ord]
    df |> mutate(
      cat1 = factor(cat1, levels = lvls),
      cat2 = factor(cat2, levels = lvls)
    )
  })() |>
  ggplot(aes(x = cat1, y = cat2, fill = med_cor)) +
  geom_tile() +
  theme_custom() +
  scale_fill_gradient2(
    high  = "#C33764",
    mid   = "white",
    low   = "#1D2671",
    midpoint = 0
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3)) +
  labs(x = "", 
       y = "",
       fill=expression("Median "~rho))

## Within Category Correlation ------------------------
cor_sum_df |>
    filter(corr_type == "Within-category") |>
    ggplot(aes(x = reorder(cat1,correlation), 
               y = correlation,
               fill = cat1)) +
    geom_boxplot(alpha = 0.7,linewidth=0.005) +
    geom_jitter(aes(color=cat1),alpha=.1)+
    theme_classic(base_size = 12) +
    scale_fill_manual(values=cat_colors)+
    scale_color_manual(values=cat_colors)+
    labs(title = "",
         x = "",
         y = "Correlation") +
    guides(color="none",
           fill="none")+geom_hline(yintercept = 0,linetype="dashed",color="grey25")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


cor_sum_df |>
    filter(corr_type == "Within-category") |>
    ggplot(aes(y = reorder(cat1,correlation), 
               x = correlation,
               fill = cat1)) +
    geom_boxplot(alpha = 0.7,linewidth=0.005) +
    geom_jitter(aes(color=cat1),alpha=.1)+
    theme_classic(base_size = 12) +
    scale_fill_manual(values=cat_colors)+
    scale_color_manual(values=cat_colors)+
    labs(title = "",
         y = "",
         x = "Correlation") +
    guides(color="none",
           fill="none")+
  geom_vline(xintercept = 0,
             linetype="dashed",
             color="grey25")+
    theme(axis.text.y = element_text(angle = 45, hjust = 1))


## Correlation with FEV1/FVC ----------------
extract_results(expom,
                result = "correlation") |>
    pluck("exposures") |> 
    filter(var1=="pftfev1fvc_actual",
           !grepl("^pft",var2)) |> 
    ggplot(aes(x = correlation, 
               y = reorder(var2,correlation),
               fill = correlation)) +
    geom_bar(stat = "identity",color="black",linewidth = 0.005)+
    scale_fill_gradient2(
      high  = "#C33764",
      mid   = "white",
      low   = "#1D2671",
        midpoint = 0
    ) +
    labs(
        x = "Correlation",
        y = NULL,
        fill = "Correlation"
    ) +
    theme_classic() +
    theme(
        legend.position = "right"
    )
      


## Correlation with PM10 ----------------
extract_results(expom,result = "correlation") |>
    pluck("exposures") |> 
    filter(var1=="pm10f",
           !grepl("^pft",var2)) |> 
    ggplot(aes(x = correlation, 
               y = reorder(var2,correlation),
               fill = correlation)) +
    geom_bar(stat = "identity",color="black",linewidth = 0.005)+
    scale_fill_gradient2(
      high  = "#C33764",
      mid   = "white",
      low   = "#1D2671",
        midpoint = 0
    ) +
    labs(
        x = "Correlation",
        y = NULL,
        fill = "Correlation"
    ) +
    theme_classic() +
    theme(
        legend.position = "right"
    )


## Sample Summary Stats ----------------------

expom |> 
    tidyexposomics::pivot_sample() |>
  mutate(gli_sex=case_when(
    gli_sex=="M" ~ "Male",
    gli_sex=="F" ~ "Female"
  ),
  black = case_when(
    black == "black" ~ "Black",
    black == "non-black" ~ "Non-Black"
  ),
  income5 = case_when(
    income5 == "RA/DK" ~ "No Response",
    .default = income5
  )) |> 
    transmute(
        Sex     = factor(gli_sex),
        `Age (yrs)`     = as.numeric(gli_age),
        `Height (cm)`  = as.numeric(gli_height),
        Ancestry       = factor(black, 
                                levels = c("Non-Black", "Black")),
        `FEV1/FVC` = pftfev1fvc_actual,
        BMI     = factor(bmicat4),
        `Household Income` = income5
    ) |> tbl_summary(
    statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    missing = "no"
) |> 
    modify_header(label = "**Characteristic**") |> 
    bold_labels()




## Association with FEV1/FVC --------------------------
lm(
    `FEV1/FVC` ~ Sex + Ancestry + BMI + `Household Income` + `Age (yrs)` + `Height (cm)`,
    data = expom |> 
        tidyexposomics::pivot_sample() |>
        mutate(gli_sex=case_when(
            gli_sex=="M" ~ "Male",
            gli_sex=="F" ~ "Female"
        ),
        black = case_when(
            black == "black" ~ "Black",
            black == "non-black" ~ "Non-Black"
        ),
        income5 = case_when(
            income5 == "RA/DK" ~ "No Response",
            .default = income5
        )) |> 
        transmute(
            Sex     = factor(gli_sex),
            `Age (yrs)`     = as.numeric(gli_age),
            `Height (cm)`  = as.numeric(gli_height),
            Ancestry       = factor(black, 
                                    levels = c("Non-Black", "Black")),
            `FEV1/FVC` = pftfev1fvc_actual,
            BMI     = factor(bmicat4),
            `Household Income` = income5
        )
) |> 
    gtsummary::tbl_regression() |> 
  bold_labels() 


## Association with Cell Types --------------------
# Perform ExWAS Analysis
cell_res <- expom |>
    run_association(
        source = "exposures",
        outcome = "pftfev1fvc_actual",
        feature_set = aw_cb |>
          filter(category=="CBC Differential") |> 
          pull(variable),
        covariates = c(
          "gli_age",
          "gli_sex",
          "gli_height",
          "fis"
        ),
        action = "get",
        family = "gaussian"
    ) |> 
  pluck("results_df") |> 
  filter(term %in% c(
    aw_cb |>
          filter(category=="CBC Differential") |> 
          pull(variable)
  ))

cell_res |> 
  # Order categories by median effect size
    mutate(
        sig = p.value< 0.05,
        term = case_when(
          sig == TRUE ~ paste0("* ",definition),
          .default =definition),
        # Order terms by estimate (within categories)
        term = forcats::fct_reorder(term, estimate)
    ) |> 
    ggplot(aes(x = estimate, 
               y = term,
               color = estimate,
               fill = estimate)) +
    geom_point(shape = 21, 
               size = 2.5, 
               stroke = 0.8) +
    geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                       xmax = estimate + 1.96 * std.error),
                   height = 0.15) +
    scale_color_gradient2(
        low="blue4",
        mid="grey85",
        high="red4",
        midpoint = 0
    )+
    scale_fill_gradient2(
        low="blue4",
        mid="grey85",
        high="red4",
        midpoint = 0
    )+
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
        title = NULL,
        x = "Effect size",
        y = "",
        fill ="Effect size",
        color ="Effect size"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle=0,
                                      face="bold.italic"))


## Cytokine association with FEV1/FVC -----------------------
# Perform ExWAS Analysis
cell_cyto_res <- expom |>
    run_association(
        source = "exposures",
        outcome = "pftfev1fvc_actual",
        feature_set = aw_cb |>
    filter(grepl(
      "Serum Cytokines|Nasal Cytokines|CBC Differential",category)) |> 
    filter(variable %in% colnames(
      tidyexposomics::pivot_sample(expom))) |> 
    pull(variable),
        covariates = c(
          "gli_age",
          "gli_sex",
          "gli_height",
          "fis"
        ),
        action = "get",
        family = "gaussian"
    ) |> 
  pluck("results_df") |> 
  filter(term %in% c(
    aw_cb |>
    filter(grepl(
      "Serum Cytokines|Nasal Cytokines|CBC Differential",
      category)) |> 
    filter(variable %in% colnames(
      tidyexposomics::pivot_sample(expom))) |> 
    pull(variable)
  ))

cell_cyto_res |> 
  # Order categories by median effect size
    mutate(
        sig = p.value< 0.05,
        term = case_when(
          sig == TRUE ~ paste0("* ",term),
          .default =term),
        # Order terms by estimate (within categories)
        term = forcats::fct_reorder(term, estimate),
        category = forcats::fct_reorder(
            category, estimate, .fun = median, .desc = TRUE),
    ) |> 
    ggplot(aes(x = estimate, 
               y = term,
               color = estimate,
               fill = estimate)) +
    geom_point(shape = 21, 
               size = 2.5, 
               stroke = 0.8) +
    geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                       xmax = estimate + 1.96 * std.error),
                   height = 0.15) +
  ggh4x::facet_grid2(category ~ ., 
                       scales = "free",
                       space = "free") +
    scale_color_gradient2(
        low="blue4",
        mid="grey85",
        high="red4",
        midpoint = 0
    )+
    scale_fill_gradient2(
        low="blue4",
        mid="grey85",
        high="red4",
        midpoint = 0
    )+
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
        title = NULL,
        x = "Effect size",
        y = "",
        fill ="Effect size",
        color ="Effect size"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle=0,
                                      face="bold.italic"))


## Which Covariates to Include? -----------------
source("./scripts/covar_screen.R")
covars <- covar_screen(
  expom |> 
    tidyexposomics::pivot_sample(),
  outcome = "pftfev1fvc_actual",
  base_covars = c(
    "gli_age",
    "gli_sex"),
  candidate_covars = c(
    "fis",
    "black",
    "gli_height",
    "bmizscore",
    "somecol",
    "income5",
    "season",
    "med_ics")
  # candidate_covars=vars[!vars %in% c(exp_vars,"gli_age","gli_sex")]
  
)

covars$recommended_covariates


## ExWAS ---------


# Perform ExWAS Analysis
fev1_fvc_expom <- expom |>
    run_association(
        source = "exposures",
        outcome = "pftfev1fvc_actual",
        feature_set = exp_vars,
        covariates = c(
          "gli_age",
          "gli_sex",
          "fis"
        ),
        log_trans = F,
        action = "add",
        family = "gaussian"
    )

x <- expom |>
  run_association(
    source = "exposures",
    outcome = "pftfev1fvc_actual",
    feature_set = exp_vars,
    covariates = c(
      "gli_age",
      "gli_sex",
      "fis",
      "med2c"
    ),
    log_trans = F,
    action = "add",
    family = "gaussian"
  ) |> 
  extract_results(result="association")
fev1_fvc_expom |>
    plot_association(
        subtitle = paste(
          "Covariates: ",
          "Age, Sex, Height, Food Insecurity"
        ,sep = "\n"),
        source = "exposures",
        terms = exp_vars,
        filter_thresh = 0.05,
        filter_col = "p.value",
        r2_col = "adj_r2"
    )+
  labs(
    caption="P < 0.05",
    title=NULL
  )+
  theme(plot.caption = element_text(face="italic"))


## Distribution of Estimates per Category --------------
assoc_df <- fev1_fvc_expom |> 
  extract_results(result = "association") |> 
  pluck("assoc_exposures",
        "results_df") |> 
  # Order categories by median effect size
  mutate(
    sig = p.value< 0.05,
    clean_name = case_when(
      sig == TRUE ~ paste0("* ",clean_name),
      .default =clean_name),
    # Order terms by estimate (within categories)
    clean_name = forcats::fct_reorder(clean_name, estimate),
    category = forcats::fct_reorder(
      category, estimate, .fun = median, .desc = TRUE),
  ) 

assoc_df |> 
ggplot(aes(
  x=reorder(category,estimate),
  y=estimate,
  fill=category,
  color=category
))+
  geom_jitter(alpha=.1)+
  geom_boxplot(color="black")+
  theme_custom()+
  theme(axis.text.x=element_text(angle=70))+
  geom_hline(yintercept = 0,
             linetype="dashed")+
  guides(color="none",
         fill="none")+
  scale_color_manual(values=cat_colors)+
  scale_fill_manual(values=cat_colors)+
  labs(
    x="",
    y="Estimate"
  )+
  coord_flip()+
  theme(axis.text.x=element_text(angle = 0,hjust=.3))

### Estimate Tile Plot ----------------
assoc_df |> 
  ggplot(aes(
    x=1,
    y=
  ))

## Manahattan Plot ------------
plot_manhattan(fev1_fvc_expom,
               facet_angle = 0,
               facet_cols = cat_colors)


## Grab the ExWAS Results and Plot ---------------

assoc_df |> 
    ggplot(aes(x = estimate, 
               y = clean_name,
               color = estimate,
               fill = estimate)) +
    geom_point(shape = 21, 
               size = 2.5, 
               stroke = 0.8) +
    geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                       xmax = estimate + 1.96 * std.error),
                   height = 0.15) +
    ggh4x::facet_grid2(category ~ ., 
                       scales = "free",
                       space = "free", 
                       strip = ggh4x::strip_themed(
                           background_y = ggh4x::elem_list_rect(
                               fill = cat_colors[levels(assoc_df$category)] |> 
                                 scales::alpha(.3)
                               ))) +
    scale_color_gradient2(
       high = "#9B6981FF",
        mid = "grey85",
        low = "#5399b0",
        midpoint = 0,
       limits = c(-.03, 0.03),
       breaks = c(-.03,0,.03)
    )+
    scale_fill_gradient2(
      high = "#9B6981FF",
      mid = "grey85",
      low = "#5399b0",
      midpoint = 0,
      limits = c(-.03, 0.03),
      breaks = c(-.03,0,.03)
    )+
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
        title = NULL,
        x = "Effect size",
        y = "",
        fill ="Effect size",
        color ="Effect size"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          #legend.direction = "vertical",
          legend.title.position = "top",
          strip.text.y = element_text(angle=0,
                                      face="bold.italic"))

## How many are significant ---------------

assoc_df |> 
  mutate(sig = p.value<0.05) |> 
  dplyr::count(sig) |> 
  mutate(sig_label = ifelse(
    sig==TRUE,
    "Significant",
    "Not Significant")) |> 
  ggplot(aes(
    x=n,
    y=reorder(sig_label,n),
    fill=sig_label
  ))+
  geom_col()+
  theme_custom()+
  scale_fill_manual(
    values = c(
      "Significant" = "#AD5A6BFF",
      "Not Significant" ="grey90"
    )
  )+
  theme(legend.position = "none")+
  labs(
    x="No. of Associations",
    y=""
  )+
  coord_flip()


## ExWAS Stability Analysis ---------------------
source("./scripts/run_exwas_sensitivity.R")
exwas_stability_res <- run_exwas_sensitivity(
    df = tidyexposomics::pivot_sample(fev1_fvc_expom),
    exposures = exp_vars,
    outcome = "pftfev1fvc_actual",
    covariates = c("gli_age" , "gli_sex", "fis"),
    covariates_to_remove = c("gli_age" , "gli_sex", "fis"),
    family = gaussian(),
    bootstrap_n = 100,
    parallel = TRUE
)

sens_df <- exwas_stability_res$sensitivity_df |>
    mutate(
        signif = (p_value <= 0.05) ,
        sign = sign(estimate)
    ) %>%
    group_by(exposure) %>%
    summarise(
        n_tests = n(),
        n_signif = sum(signif, na.rm = TRUE),
        prop_signif = n_signif / n_tests,
        median_effect = median(estimate, na.rm = TRUE),
        median_p=median(p_value, na.rm = TRUE),
        sd_effect = sd(estimate, na.rm = TRUE),
        sign_consistency = ifelse(
          all(is.na(sign)), NA_real_,
          max(mean(sign > 0, na.rm = TRUE), 
              mean(sign < 0, na.rm = TRUE))
        ),
        stability_score = prop_signif * sign_consistency,
        .groups = "drop"
    ) |> 
  left_join(
    aw_cb,
    by=c("exposure"="variable")
  )




## Make the Sensitivity Analysis Summary Df --------------
sens_assoc_df <- sens_df |> 
    # Order categories by median effect size
    mutate(
        category = forcats::fct_reorder(
            category, median_effect, .fun = median, .desc = TRUE),
        
        # Order terms by estimate (within categories)
        exposure = forcats::fct_reorder(exposure, median_effect),
        
        # Significance flag for fill aesthetics
        sig = median_p < 0.05
    ) 
sens_assoc_df |> 
    ggplot(aes(x = median_effect, 
               y = exposure,
               color = median_effect,
               fill = median_effect)) +
    geom_point(shape = 21, 
               size = 2.5, 
               stroke = 0.8) +
    geom_errorbarh(aes(xmin = median_effect - 1.96 * sd_effect,
                       xmax = median_effect + 1.96 * sd_effect),
                   height = 0.15) +
    #facet_grid(category~., scales = "free",space = "free") +
    ggh4x::facet_grid2(category ~ ., 
                       scales = "free",
                       space = "free", 
                       strip = ggh4x::strip_themed(
                           background_y = ggh4x::elem_list_rect(
                               fill = cat_colors[levels(sens_assoc_df$category)]))) +
    # ggh4x::force_panelsizes(rows = c(
    #   2,2,2,2,2,1,1,2,1,1,3,3,3
    # )) + 
    scale_color_gradient2(
        low="blue4",
        mid="grey85",
        high="red4",
        midpoint = 0
    )+
    scale_fill_gradient2(
        low="blue4",
        mid="grey85",
        high="red4",
        midpoint = 0
    )+
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
        title = NULL,
        x = "Effect size",
        y = "",
        fill ="Effect size",
        color ="Effect size"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle=0,
                                      face="bold.italic"))


## Sign Consistency ---------------------
sign_cons_df <- sens_assoc_df |> 
  mutate(
    clean_name = tidytext::reorder_within(x = clean_name,by = sign_consistency,within = category
  ))  
  



sign_cons_df <- sign_cons_df |> 
  mutate(category = factor(category, levels = c(
    sign_cons_df |> 
      group_by(category) |> 
      reframe(cat_cons = mean(sign_consistency)) |> 
      arrange(desc(cat_cons)) |> 
      pull(category)
  )))

sign_cons_df |>
  ggplot(aes(x = sign_consistency, 
             y = clean_name,
             color = sign_consistency,
             fill = sign_consistency)) +
  geom_point(shape = 21, 
             size = 2.5, 
             stroke = 0.8) +
  geom_segment(aes(
    x = 0,
    xend = sign_consistency,
    y=clean_name
  ))+
  tidytext::scale_y_reordered()+
  ggh4x::facet_grid2(category ~ ., 
                     scales = "free",
                     space = "free", 
                     strip = ggh4x::strip_themed(
                       background_y = ggh4x::elem_list_rect(
                         fill = cat_colors[levels(sign_cons_df$category)] |> 
                           scales::alpha(.3)
                         
                         ))) +
  # ggh4x::force_panelsizes(rows = c(
  #   2,2,2,2,2,1,1,2,1,1,3,3,3
  # )) + 
  scale_color_gradientn(
    colors = c("#000004FF","#B63679FF","#FCFDBFFF") |> rev()
  )+
  scale_fill_gradientn(
    colors = c("#000004FF","#B63679FF","#FCFDBFFF") |> rev()
  )+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = NULL,
    x = "Sign Consistency",
    y = "",
    fill ="Sign Consistency",
    color ="Sign Consistency"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title.position = "top",
        strip.text.y = element_text(angle=0,
                                    face="bold.italic"))


## How many are over 75% consistent? -------------

sens_assoc_df |> 
  mutate(over_75_cons = ifelse(sign_consistency > .75,
                                "Consistent",
                                "Not Consistent")) |> 
  dplyr::count(over_75_cons) |> 
  ggplot(aes(
    x=n,
    y=reorder(over_75_cons,n),
    fill=over_75_cons
  ))+
  geom_col()+
  theme_custom()+
  scale_fill_manual(
    values = c(
      "Consistent" = "#19547b",
      "Not Consistent" ="grey90"
    )
  )+
  theme(legend.position = "none")+
  labs(
    x="No. of Associations",
    y=""
  )+
  coord_flip()


## Median Estimate -------------------------
sens_assoc_df |> 
    group_by(category) |> 
    mutate(med_cat=mode(median_effect)) |> 
    ungroup() |> 
    ggplot(aes(
        x = reorder(category,med_cat),
        y = median_effect,
        fill = category,
        color = category)) +
    geom_boxplot(outlier.shape = NA, 
                 alpha = 0.8,
                 color = "black",
                 linewidth=0.005) +
    geom_jitter(alpha=0.1)+
    labs(
        y = "Median Estimate",
        x = "",
        title = ""
    ) +
    scale_fill_manual(values=cat_colors)+
    scale_color_manual(values=cat_colors)+
    theme_classic(base_size = 12) +
    theme(legend.position = "none")+
    rotate_x_text(angle = 65)+
    ylim(c(-0.04,0.04))+
    geom_hline(yintercept = 0,
               linetype="dashed")


## Median Sign Consistency -------------------------
sens_assoc_df |> 
    group_by(category) |> 
    mutate(med_cat=median(sign_consistency)) |> 
    ungroup() |> 
    ggplot(aes(
        x = reorder(category,med_cat),
        y = sign_consistency,
        fill = category,
        color = category)) +
    geom_boxplot(outlier.shape = NA, 
                 alpha = 0.8,
                 color = "black") +
    geom_jitter(alpha=0.1)+
    labs(
        y = "Median Sign Consistency",
        x = "",
        title = ""
    ) +
    scale_fill_manual(values=cat_colors)+
    scale_color_manual(values=cat_colors)+
    theme_classic(base_size = 12) +
    theme(legend.position = "none")+
    rotate_x_text(angle = 65)
    # geom_hline(yintercept = 0.5,
    #            linetype="dashed")


## Power Analysis ----------------------
source("./scripts/sample_size.R")

# grab association results
assoc_df <- fev1_fvc_expom@metadata |> 
  pluck(
    "association",
    "assoc_exposures",
    "results_df"
  )

# n = your sample size, covariates = "gli_age" + "gli_sex" + "gli_height" +"fis" = 4
sample_size_results <- bind_rows(
  calc_adj_sample_sizes(assoc_df, 
                        model_n = 48, 
                        n_covariates = 4, 
                        alpha = 0.05, 
                        power = 0.70) |>
    mutate(power = "0.70"),
  calc_adj_sample_sizes(assoc_df,
                        model_n = 48,
                        n_covariates = 4,
                        alpha = 0.05, 
                        power = 0.80) |>
    mutate(power = "0.80"),
  calc_adj_sample_sizes(assoc_df, 
                        model_n = 48,
                        n_covariates = 4,
                        alpha = 0.05, 
                        power = 0.90) |>
    mutate(power = "0.90")
)

sample_size_results |> 
  ggplot(aes(
    x = required_n,
    y = abs(estimate),
    fill = power
  ))+
  geom_point(alpha=0.7,size=2.5,shape=21,color="black")+
  geom_line(alpha=0.5)+
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  theme_custom()+
  scale_fill_manual(values=c(
    "#483d5b",
    "#855174",
    "#c3687a"
  ))+
  labs(
    x="Sample Size",
    y="|Estimate|",
    fill="Power"
  )

## GAM for Non-Linearity ---------------
source("./scripts/run_gam.R")
gam_res <- run_gam(
  df = pivot_sample(fev1_fvc_expom),
  codebook = fev1_fvc_expom |>
    extract_results(result = "codebook"),
  outcome = "pftfev1fvc_actual",
  feature_set = exp_vars,
  covariates = c(
    "gli_age",
    "gli_sex",
    "gli_height",
    "fis"
  )
)

## How many are non-linear? ---------------
gam_res |> 
  mutate(non_lin = ifelse(edf>2 & p.value<0.05,
                          "Non-Linear",
                          "Linear")) |> 
  dplyr::count(non_lin) |> 
  ggplot(aes(
    x=n,
    y=reorder(non_lin,n),
    fill=non_lin
  ))+
  geom_col()+
  theme_custom()+
  scale_fill_manual(
    values = c(
      "Non-Linear" = "#360033",
      "Linear" ="grey90"
    )
  )+
  theme(legend.position = "none")+
  labs(
    x="No. of Associations",
    y=""
  )+
  coord_flip()

## Distribution EDF by Category --------------
gam_res |> 
  ggplot(aes(
    x=reorder(category,edf),
    y=edf,
    fill=category,
    color=category
  ))+
  geom_jitter(alpha=.1)+
  geom_boxplot(color="black")+
  theme_custom()+
  theme(axis.text.x=element_text(angle=70))+
  geom_hline(yintercept = 1,
             linetype="dashed")+
  guides(color="none",
         fill="none")+
  scale_color_manual(values=cat_colors)+
  scale_fill_manual(values=cat_colors)+
  labs(
    x="",
    y="EDF"
  )+
  coord_flip()+
  theme(axis.text.x=element_text(angle = 0,hjust=.3))
## Elastic Net ------------------
set.seed(42)
enet_res <-
  fev1_fvc_expom |>
  tidyexposomics::pivot_sample() |>
  select(pftfev1fvc_actual, 
         all_of(exp_vars), 
         all_of(c("gli_age","gli_sex","gli_height","fis"))) |>
  na.omit() |>
  (\(df) {
    
    X <- model.matrix(
      ~ . - pftfev1fvc_actual - 1,
      data = df
    )
    
    y <- df$pftfev1fvc_actual
    
    penalty <- ifelse(colnames(X) %in% exp_vars, 1, 0)
    
    cv.glmnet(
      x = X,
      y = y,
      family = "gaussian",
      alpha = 0.5,                 # elastic net
      penalty.factor = penalty,
      standardize = TRUE,
      nfolds = 10
    )
  })()

enet_tbl <-
  broom::tidy(enet_res$glmnet.fit) |>
  filter(lambda==enet_res$lambda.1se)


## Quantile G-Computation -----

library(qgcomp)


gcomp_res <- map(fev1_fvc_expom@metadata$codebook |>
      filter(variable %in% num_exp_vars) |> 
      pull(category) |> 
      unique(),~{
        
        message("Working on category: ",.x)
        gcomp_vars <- fev1_fvc_expom@metadata$codebook |>
          filter(variable %in% num_exp_vars) |>
          filter(category == .x) |> 
          pull(variable)
        
        covars <- c(
          "gli_age",
          "gli_sex",
          "fis"
        )
        
        outcome <- "pftfev1fvc_actual"
        
        fit <- qgcomp.noboot(
          pftfev1fvc_actual ~ . + gli_age + gli_sex + fis,
          expnms  = gcomp_vars,
          data    = pivot_sample(fev1_fvc_expom) |> 
            dplyr::select(
              all_of(covars),
              all_of(outcome),
              all_of(gcomp_vars)
            ),
          q       = 4,
          bayes   = TRUE,
          family  = gaussian()
        ) |> 
          broom::tidy() |> 
          filter(term != "(Intercept)") |> 
          mutate(category=.x)
      }) |> 
  bind_rows()


## Weighted Quantile Sum -----

library(gWQS)

wqs_res <- map(
  fev1_fvc_expom@metadata$codebook |>
    filter(variable %in% num_exp_vars) |> 
    pull(category) |> 
    unique(),
  ~{
    
    message("Working on category: ",.x)
    wqs_vars <- fev1_fvc_expom@metadata$codebook |>
      filter(variable %in% num_exp_vars) |>
      filter(category == .x) |> 
      pull(variable)
    
    covars <- c(
      "gli_age",
      "gli_sex",
      "fis"
    )
    
    outcome <- "pftfev1fvc_actual"
    
    fit <- gwqs(
      pftfev1fvc_actual ~ wqs + gli_age + gli_sex + fis,
      mix_name = wqs_vars,
      data     = pivot_sample(fev1_fvc_expom) |> 
        dplyr::select(
          all_of(covars),
          all_of(outcome),
          all_of(wqs_vars)
        ),
      q        = 4,   
      validation = 0.7,
      b        = 100,
      family   = "gaussian"
    ) |> 
      pluck("fit") |> 
      broom::tidy() |> 
      filter(term =="wqs") |> 
      mutate(category=.x)
  }) |> 
  bind_rows()



## Save Data ----------------------
# Save MultiAssayExperiment Objects
saveRDS(fev1_fvc_expom,file = "./results/exwas/fev1_fvc_expom.rds")

# Save Bootstrapping Results
saveRDS(exwas_stability_res,file =
          "./results/exwas/exwas_stability_res.rds")


