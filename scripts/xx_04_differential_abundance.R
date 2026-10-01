## Load Libraries ------------------------

library(tidyverse)
library(janitor)
library(ggpubr)
library(furrr)
library(ggsci)
library(tidyexposomics)
library(MultiAssayExperiment)
library(GSVA)
library(compositions)
library(broom)

source("./scripts/bin/internals.R")
source("./scripts/bin/labels_clean.R")

## Load Data --------------------
# Codebook
aw_cb <- readRDS("./results/input_data/aw_cb.rds")

# MultiAssayExperiment Object
fev1_fvc_expom <- readRDS("./results/exwas/fev1_fvc_expom.rds")

# Exposure variables
exp_vars <- readRDS("./results/qc/exp_vars.rds")
num_exp_vars <- readRDS("./results/qc/num_exp_vars.rds")

# All variables
vars <- aw_cb$variable |> 
  (\(chr)chr[chr %in% colnames(colData(fev1_fvc_expom))])()

# Setting colors for exposure categories
cat_colors <- scales::alpha(get_palette("aaas",13), 0.5) |> 
  (\(colors){
    names(colors) <- aw_cb |> 
      filter(variable %in% exp_vars) |> 
      pull(category) |> 
      unique()
    colors
  })()

msig <- msigdbr::msigdbr(species = "human",collection = "C5")

## Normalize ------------------
normalize_rna_assay <- function(mae, assay_name) {
  se <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  mat <- SummarizedExperiment::assay(se, "counts") |> as.matrix()

  dge <- edgeR::DGEList(counts = mat) |>
    edgeR::calcNormFactors(method = "TMM")

  log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

  # shift to non-negative to avoid the pseudocount patch in run_differential_abundance
  # this preserves relative differences, which is all limma_trend needs
  log_cpm_shifted <- log_cpm - min(log_cpm, na.rm = TRUE)

  SummarizedExperiment::assay(se, "counts") <- log_cpm_shifted
  MultiAssayExperiment::experiments(mae)[[assay_name]] <- se
  mae
}

rna_assays <- c("CD4 T cell RNA", "CD16 Monocyte RNA",
                "CD4 T cell Isoforms", "CD16 Monocyte Isoforms")

fev1_fvc_expom <- rna_assays |>
  purrr::reduce(normalize_rna_assay, .init = fev1_fvc_expom)

## Checking Counts -----------------
fev1_fvc_expom |>
  MultiAssayExperiment::experiments() |>
  purrr::map(
    \(se) SummarizedExperiment::assay(se, "counts") |> 
      as.matrix() |>
      max(na.rm = TRUE)
  )

## Run Diff. Abd. ------------------
source("./scripts/run_sva_limma_trend.R")
# Run differential abundance analysis
# run across RNA assays only, bind results
base_formula <- ~ pftfev1fvc_actual + gli_age + gli_sex + fis

rna_assays <- c("CD4 T cell RNA", "CD16 Monocyte RNA",
                "CD4 T cell Isoforms", "CD16 Monocyte Isoforms")

da_res_sva <- rna_assays |>
  purrr::map(run_sva_limma_trend, 
             mae = fev1_fvc_expom,
             base_formula = base_formula) |>
  purrr::list_rbind()

da_res_reg <- map(c(
  "Proteomics","Adductomics","Adduct Burden" 
),~{
  se  <- tidyexposomics:::.update_assay_colData(fev1_fvc_expom, .x)
  # --- limma-trend ---
  tidyexposomics:::.run_limma_trend(
    se             = se,
    formula        = base_formula,
    abundance_col  = "counts",
    scaling_method = "none"
  ) |>
    dplyr::mutate(exp_name = .x)
  
}) |>
  purrr::list_rbind()

da_res <- bind_rows(
  da_res_sva,
  da_res_reg
)


# checking
da_res |> 
  group_by(exp_name) |>
  summarise(
    median_logFC = median(logFC),
    mean_logFC = mean(logFC)
  ) |>
  arrange(desc(median_logFC))

## Checking Direction ---------------
da_res |> 
  filter(P.Value<0.05) |> 
  filter(feature == "ENSG00000063587",
         exp_name == "CD16 Monocyte RNA") |> 
  pull(logFC)

fev1_fvc_expom |>  
  pivot_exp(exp_name ="CD16 Monocyte RNA",
            features ="ENSG00000063587") |>
  ggplot(aes(x=log2(counts+1),
             y=pftfev1fvc_actual))+
  geom_point()+
  theme_custom()+
  stat_smooth(method=lm)


da_res |> 
  filter(P.Value<0.05) |> 
  filter(feature == "ENSG00000214455",
         exp_name == "CD4 T cell RNA") |> 
  pull(logFC)

fev1_fvc_expom |>  
  pivot_exp(exp_name ="CD4 T cell RNA",
            features ="ENSG00000214455") |>
  ggplot(aes(x=log2(counts+1),
             y=pftfev1fvc_actual))+
  geom_point()+
  theme_custom()+
  stat_smooth(method=lm)

## No. of DEGs ----------------

da_res |>
  filter(P.Value < 0.05) |> 
  dplyr::count(exp_name,sign(logFC)) |>
  mutate(direction=ifelse(
    `sign(logFC)`>0,
    "Upregulated",
    "Downregulated")) |> 
  group_by(exp_name) |> 
  mutate(total=sum(n)) |> 
  ungroup() |> 
  ggplot(aes(
    x=n,
    y=reorder(exp_name,total),
    fill=direction
  ))+
  geom_col(alpha = 0.7,
           colour = "black") +
  scale_fill_manual(values = c(
    "#5399b0","#9B6981FF"
  ))+
  theme_custom()+
  theme(
    plot.title = element_text(face = "bold.italic"),
    plot.subtitle = element_text(face = "italic"),
    plot.caption = element_text(face = "italic")
  ) +
  scale_y_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD4 T cell Isoforms"  = expression(CD4^"+" ~ "T cell Isoforms"),
    "CD4 T cell miRNA"     = expression(CD4^"+" ~ "T cell miRNA"),
    
    "CD16 Monocyte RNA"       = expression(CD16^"+" ~ "Monocyte RNA"),
    "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
    
    "Adductomics"       = "Adductomics",
    "Adductomics Load"  = "Adductomics Load",
    "Proteomics"        = "Proteomics"
  ))+
  labs(
    x = "No. of Differential Features",
    y = "",
    fill = "Direction",
    title = "FEV1/FVC",
    subtitle = "",
    caption = "Nominal P < 0.05"
  )


## What is the relationship between omics? ----------------

pc1_res <- map(
  names(experiments(fev1_fvc_expom)),~{
    fev1_fvc_expom[[.x]] |> 
      assay() |> 
      (\(mat) log2(mat+1))() |> 
      t() |> 
      prcomp() |> 
      (\(pr) pr[["x"]])() |> 
      as.data.frame() |> 
      rownames_to_column(".sample") |> 
      dplyr::select(.sample,PC1) |> 
      mutate(exp_name =.x)
  }
) |> 
  bind_rows()

pc1_res |> 
  pivot_wider(names_from = exp_name, 
              values_from = PC1) |> 
  column_to_rownames(".sample") |> 
  cor(method = "spearman",
      use = "complete.obs") |> 
  (\(mat) {
    ord <- mat |> 
      pheatmap::pheatmap(silent = TRUE) |> _$tree_row$order
    mat |> 
      as.data.frame() |> 
      rownames_to_column("sample1") |> 
      pivot_longer(-sample1,
                   names_to = "sample2",
                   values_to = "r") |> 
      mutate(across(c(sample1, sample2),
                    \(x) factor(x, levels = rownames(mat)[ord])))
  })() |> 
  ggplot(aes(x = sample1, y = sample2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "#1D2671",
    mid      = "white",
    high     = "#C33764",
    midpoint = 0,
    limits   = c(-1, 1),
    breaks   = c(-1, -0.5, 0, 0.5, 1),
    name     = expression(~rho)
  ) +
  theme_custom() +
  scale_y_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD4 T cell Isoforms"  = expression(CD4^"+" ~ "T cell Isoforms"),
    "CD4 T cell miRNA"     = expression(CD4^"+" ~ "T cell miRNA"),
    
    "CD16 Monocyte RNA"       = expression(CD16^"+" ~ "Monocyte RNA"),
    "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
    
    "Adductomics"       = "Adductomics",
    "Adductomics Load"  = "Adductomics Load",
    "Proteomics"        = "Proteomics"
  ))+
  scale_x_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD4 T cell Isoforms"  = expression(CD4^"+" ~ "T cell Isoforms"),
    "CD4 T cell miRNA"     = expression(CD4^"+" ~ "T cell miRNA"),
    
    "CD16 Monocyte RNA"       = expression(CD16^"+" ~ "Monocyte RNA"),
    "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
    
    "Adductomics"       = "Adductomics",
    "Adductomics Load"  = "Adductomics Load",
    "Proteomics"        = "Proteomics"
  ))+
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.3),
    axis.title   = element_blank()
  )
## Overlaps -------------------
source("./scripts/get_pairwise_overlaps.R")

da_res <- x |> 
  extract_results("differential_analysis") |> 
  pluck("differential_abundance")  

da_res |> 
  filter(P.Value<0.05) |> 
  (\(df)split(df,df$exp_name))() |> 
  map(~.x |> 
        pull(feature_map) |> 
        unique()) |> 
  get_pairwise_overlaps() |> 
  arrange(jaccard) |> 
  tidygraph::as_tbl_graph() |> 
  ggraph::ggraph(
        layout = "linear",
        circular = TRUE
    ) +
        ggraph::geom_edge_arc(aes(
            color = jaccard,
            width = log10(num_shared+1)*2
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
        high = "#0f0c29",
        mid = "white",
        low = "#43C6AC",
        midpoint = 0.1,
        guide = ggraph::guide_edge_colorbar()
    )+
        ggplot2::coord_fixed(xlim = c(-2, 2), ylim = c(-2, 2)) +
        ggplot2::guides(edge_width = "none") +
        ggplot2::labs(
            edge_color = "Jaccard Coef."
        )


## LogFC ----------------
fev1_fvc_expom |> 
  extract_results("differential_analysis") |> 
  pluck("differential_abundance") |> 
  filter(P.Value<0.05) |> 
  group_by(exp_name) |> 
  mutate(mean_abs_logfc=mean(abs(logFC))) |> 
  ungroup() |> 
  ggplot(aes(
    x=reorder(exp_name,mean_abs_logfc),
    y=abs(logFC),
    fill=exp_name,
    color=exp_name
  ))+
  geom_boxplot(alpha=.5) +
  theme_pubr(legend = "none")+
  rotate_x_text(angle=45)+
  ggsci::scale_fill_npg()+
  ggsci::scale_color_npg()+
  labs(
    x="",
    y=expression("Abs(Log"[2]*"FC)")
  )


## Sensitivity Analysis ---------------------
# Perform Sensitivity Analysis
fev1_fvc_expom <- fev1_fvc_expom |>
    run_sensitivity_analysis(
        base_formula = ~ pftfev1fvc_actual + gli_age + gli_sex + gli_height + fis,
        methods = c("limma_trend"),
        scaling_methods = c("none"),
        covariates_to_remove = c(
            "gli_age",
            "gli_sex",
            "gli_height",
            "fis"
        ),
        pval_col = "P.Value",
        logfc_col = "logFC",
        logFC_threshold = log2(1),
        pval_threshold = 0.05,
        stability_metric = "stability_score",
        bootstrap_n = 100,
        action = "add"
    )


x <- fev1_fvc_expom |>
  run_sensitivity_analysis_fast(
    base_formula = ~ pftfev1fvc_actual + gli_age + gli_sex + gli_height + fis,
    methods = c("limma_trend"),
    scaling_methods = c("none"),
    covariates_to_remove = c(
      "gli_age",
      "gli_sex",
      "gli_height",
      "fis"
    ),
    da_res = da_res |> filter(P.Value<0.05),
    pval_col = "P.Value",
    logfc_col = "logFC",
    pval_threshold = 0.05,
    bootstrap_n = 500,
    action = "add"
  )

# fev1_fvc_expom@metadata$differential_analysis$sensitivity_analysis$sensitivity_df <- NULL

## No. of DEGs (bootstrapped) ----------------

da_res <- x |> 
  extract_results("differential_analysis") |> 
  pluck("differential_abundance") 

x |> 
  extract_results("differential_analysis") |> 
  pluck(
    "sensitivity_analysis",
    "feature_stability"
  ) |> 
  inner_join(da_res,
             by=c("feature","exp_name")) |>
  filter(sign_consistency>.75,prop_sig>.5) |>
  dplyr::count(exp_name,sign(logFC)) |>
  mutate(direction=ifelse(
    `sign(logFC)`>0,
    "Upregulated",
    "Downregulated")) |> 
  group_by(exp_name) |> 
  mutate(total=sum(n)) |> 
  ungroup() |> 
  ggplot(aes(
    x=n,
    y=reorder(exp_name,total),
    fill=direction
  ))+
  geom_col(alpha = 0.7,
           colour = "black") +
  scale_fill_manual(values = c(
    "#5399b0","#9B6981FF"
  ))+
  theme_custom()+
  theme(
    plot.title = element_text(face = "bold.italic"),
    plot.subtitle = element_text(face = "italic"),
    plot.caption = element_text(face = "italic")
  ) +
  scale_y_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD4 T cell Isoforms"  = expression(CD4^"+" ~ "T cell Isoforms"),
    "CD4 T cell miRNA"     = expression(CD4^"+" ~ "T cell miRNA"),
    
    "CD16 Monocyte RNA"       = expression(CD16^"+" ~ "Monocyte RNA"),
    "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
    
    "Adductomics"       = "Adductomics",
    "Adductomics Load"  = "Adductomics Load",
    "Proteomics"        = "Proteomics"
  ))+
  labs(
    x = "No. of Differential Features",
    y = "",
    fill = "Direction",
    title = "FEV1/FVC",
    subtitle = "",
    caption = "Nominal P < 0.05"
  )

## Go Enrichment --------------
source("./scripts/gsea_go_enrichment.R")
# Get GO Data
go_term_data <- {
  options(GO_ANNOTATION_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN")
  
  go_mapping <- fenr:::fetch_go_genes_go(
    species   = "goa_human",
    use_cache = TRUE,
    on_error  = "stop"
  )
  
  go_terms <- fenr:::fetch_go_terms(
    use_cache = TRUE,
    on_error  = "stop"
  )
  
  list(terms = go_terms, mapping = go_mapping)
}

da_lst <- da_res |> 
  filter(P.Value<0.05) |> 
  mutate(direction=ifelse(
    logFC>0,
    "Upregulated",
    "Downregulated")) |>
  mutate(combo = paste(exp_name,direction,sep="_")) |> 
  mutate(feature_map = ifelse(is.na(feature_map),
                              feature_id,
                              feature_map)) |> 
  (\(df) split(df,df$combo))() |> 
  map(~{
    .x |> 
      pull(feature_map) |> 
      unique()
  }) |> 
  (\(lst) lst[lengths(lst)>10])()

# Build universes per cell type by modality
universes <- pivot_feature(fev1_fvc_expom) |>
  filter(.exp_name %in% c(
    "CD4 T cell RNA", "CD4 T cell Isoforms",
    "CD16 Monocyte RNA", "CD16 Monocyte Isoforms"
  )) |>
  dplyr::count(.exp_name, feature_map) |>   # deduplicates
  split(~.exp_name) |>
  map(~ pull(.x, feature_map))

# Helper to pick the right universe from the group name
pick_universe <- function(group_name) {
  exp_key <- case_when(
    grepl("CD4.*RNA",      group_name) ~ "CD4 T cell RNA",
    grepl("CD4.*Isoform",  group_name) ~ "CD4 T cell Isoforms",
    grepl("CD16.*RNA",     group_name) ~ "CD16 Monocyte RNA",
    grepl("CD16.*Isoform", group_name) ~ "CD16 Monocyte Isoforms"
  )
  universes[[exp_key]]
}

# Run enrichment with matched universes
enr_res <- future_imap(
  da_lst,
  ~ {
    options(GO_ANNOTATION_URL = "http://current.geneontology.org/annotations/goa_human.gaf.gz")
    run_fenr(
      selected_genes  = .x,
      term_data       = go_term_data,
      universe_genes  = pick_universe(.y),
      db              = "GO",
      species         = "goa_human",
      feature_col     = "gene_symbol"
    )
  },
  .progress = TRUE
) |>
  bind_rows(.id = "group")

enr_res_sig <- enr_res |>
  filter(
    p_adjust < .1,
    n_with_sel > 5,
    N_with < 1000
  )

## GSEA ---------------------

source("./scripts/gsea_go_enrichment.R")

# Get GO Data
go_term_data <- {
  options(GO_ANNOTATION_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN")
  
  go_mapping <- fenr:::fetch_go_genes_go(
    species   = "goa_human",
    use_cache = TRUE,
    on_error  = "stop"
  )
  
  go_terms <- fenr:::fetch_go_terms(
    use_cache = TRUE,
    on_error  = "stop"
  )
  
  list(terms = go_terms, mapping = go_mapping)
}

ranks <- da_res |>
  filter(!grepl("Isoform|Adduct",exp_name)) |> 
  filter(!is.na(feature_id)) |> 
  mutate(rank_stat = sign(logFC) * -log10(P.Value)) |>
  filter(
    is.finite(rank_stat),
    !is.na(rank_stat)
  ) |>
  (\(df) split(df, df$exp_name))() |>
  purrr::map(\(grp) {
    grp |>
      slice_max(abs(rank_stat), 
                by = feature_id, 
                n = 1, 
                with_ties = FALSE) |>
      distinct(feature_id, rank_stat) |>
      tibble::deframe()
  }) |>
  (\(lst) lst[lengths(lst) > 10])()

gsea_res <- map(
  ranks,
  \(ranks) run_fgsea(
    ranked_genes = ranks,
    db           = "GO",
    term_data    = go_term_data,
    species      = "goa_human",
    feature_col  = "gene_symbol"
  )
) |>
  dplyr::bind_rows(.id = "group") |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", ")) |>
  filter(
    size > 5,
    padj < 0.1
  )
### Plotting --------------

gsea_res |> 
  filter(size<100) |> 
  group_by(group) |> 
  arrange(desc(abs(NES))) |> 
  slice_head(n=10) |> 
  ungroup() |> 
  mutate(term_name = str_trunc(term_name,width = 40)) |> 
  mutate(term_name=tidytext::reorder_within(
    term_name,
    by=abs(NES),
    within = group)) |> 
  ggplot(aes(
    x=abs(NES),
    y=term_name,
    fill=NES
  ))+
  geom_col()+
  theme_custom()+
  tidytext::scale_y_reordered()+
  facet_grid(group~.,space="free",scales="free")+
  scale_fill_gradient2(
    high = "#9B6981FF",
    mid = "grey85",
    low = "#5399b0",
    midpoint = 0
  )+
  labs(
    x="|NES|",
    fill="NES",
    y=""
  )
## Deconvolution -------------------------------
source("./scripts/my_cibersort_win.R")
monaco <- read_tsv("./data/rna_immune_cell_monaco.tsv")

monaco_cd4_ref <- monaco |> 
  filter(`Immune cell` %in% c(
    "Memory CD4 T-cell TFH",
    "Memory CD4 T-cell Th1" ,
    "Memory CD4 T-cell Th1/Th17",
    "Memory CD4 T-cell Th17",
    "Memory CD4 T-cell Th2",
    "naive CD4 T-cell",
    "T-reg",
    "Terminal effector memory CD4 T-cell"
  )) |> 
  mutate(log=log2(nTPM+1)) |> 
  dplyr::select(Gene,`Immune cell`,log) |> 
  pivot_wider(names_from = "Immune cell",
              values_from = "log") |> 
  column_to_rownames("Gene")

cs_mat <- assay(tidyexposomics:::.update_assay_colData(
  fev1_fvc_expom,
  "CD4 T cell RNA"), "counts")

common_genes <- intersect(
  rownames(cs_mat),
  rownames(monaco_cd4_ref)
)

monaco_cd4_ref <- monaco_cd4_ref[common_genes,]

# Keep only genes that are highly expressed in at least one cell type
monaco_cd4_ref_filtered <- monaco_cd4_ref |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  filter(if_any(-gene, \(x) x > quantile(x, 0.75))) |>  # top 5% expressed in any cell type
  column_to_rownames("gene")

dim(monaco_cd4_ref_filtered)  # should be a few hundred genes

cs_mat <- cs_mat[rownames(monaco_cd4_ref_filtered),]

cs_res <- my_cibersort_win(
  Y = cs_mat,
  X = monaco_cd4_ref_filtered,
  QN   = TRUE,   
  perm = 0,      
  workers = parallel::detectCores() - 1
) |> 
  pluck("proportions") |> 
  as.data.frame()

cs_res_clean <- cs_res |>
  dplyr::select(-c(`P-value`, Correlation, RMSE)) |>
  rownames_to_column(".sample") |>
  mutate(across(-`.sample`, as.numeric))



decon <- tidyexposomics:::.update_assay_colData(
  fev1_fvc_expom,
  "CD4 T cell RNA") |>
  update_coldata(df = cs_res_clean)

## CLR testing ---------------

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
    pftfev1fvc_actual,
    gli_age,
    gli_sex,
    gli_height,
    fis
  )

model_df <- bind_cols(meta, as.data.frame(clr_mat))

# Fit one model per CLR component
decon_res <- map_dfr(
  colnames(clr_mat),
  function(clr_col) {
    lm(
      as.formula(paste(
        "pftfev1fvc_actual ~",
        clr_col,
        "+ gli_age + gli_sex + fis"
      )),
      data = model_df
    ) |>
      tidy() |>
      filter(term == clr_col) |>
      mutate(cell = sub("^clr_", "", clr_col))
  }
)


decon_res <- decon_res |>
  mutate(
    cell_clean = case_when(
      cell == "Memory.CD4.T.cell.TFH"        ~ "Tfh (memory CD4)",
      cell == "Memory.CD4.T.cell.Th1"        ~ "Th1 (memory CD4)",
      cell == "Memory.CD4.T.cell.Th2"        ~ "Th2 (memory CD4)",
      cell == "Memory.CD4.T.cell.Th17"       ~ "Th17 (memory CD4)",
      cell == "Memory.CD4.T.cell.Th1.Th17"   ~ "Th1/Th17 (memory CD4)",
      cell == "Terminal.effector.memory.CD4.T.cell" ~ "Terminal effector memory CD4",
      cell == "naive.CD4.T.cell"             ~ "Naive CD4",
      cell == "T.reg"                        ~ "Treg",
      TRUE                                   ~ cell
    )
  )

### Plot Deconvolution results -------------
decon_res |> 
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

# Adduct/Proteins -----------------
da_res |> 
  filter(grepl("Adduct Burden",exp_name)) |>
  # mutate(feature_id=ifelse(exp_name == "Adductomics",
  #                          feature_map,
  #                          feature_id)) |> 
  filter(P.Value<0.05) |> 
  view()

da_res |> 
  filter(grepl("Prote",exp_name)) |>
  filter(P.Value<0.05) |> 
  view()

### Adduct Volcano plot ------------------

a_plot <- da_res_total |>
  filter(grepl("Adduct Burden", exp_name)) |>
  mutate(
    sig = P.Value < 0.05,
    direction = case_when(
      P.Value < 0.05 & logFC > 0 ~ "Up",
      P.Value < 0.05 & logFC < 0 ~ "Down",
      .default = "ns"
    ) |> factor(levels = c("Up", "ns", "Down"))
  ) |>
  ggplot(aes(
    x     = logFC,
    y     = -log10(P.Value),
    color = direction
  )) +
  geom_point(
    size =2
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0,            linetype = "dashed", color = "grey50") +
  ggrepel::geom_label_repel(
    data = \(d) d |>
      filter(direction == "Down") |>
      slice_min(logFC, n = 5),
    aes(label = feature_id),
    color              = "#5399b0",
    fill               = "white",
    size               = 3,
    fontface           = "italic",
    box.padding        = 0.5,
    min.segment.length = 0,
    seed               = 42
  ) +
  scale_color_manual(values = c(
    "Up"   = "#9B6981FF",
    "ns"   = "grey75",
    "Down" = "#5399b0"
  )) +
  labs(
    title = "Adduct Burden",
    x     = expression("Log"[2]*"FC"),
    y     = expression("-Log"[10]*"P"),
    color = "Direction"
  ) +
  theme_custom()+
  theme(
    axis.text.x=element_text(angle=0,hjust=.6)
  )

### Proteomics Volcano plot ------------------

p_plot <- da_res_total |>
  filter(grepl("Proteomics", exp_name)) |>
  mutate(
    sig = P.Value < 0.05,
    direction = case_when(
      P.Value < 0.05 & logFC > 0 ~ "Up",
      P.Value < 0.05 & logFC < 0 ~ "Down",
      .default = "ns"
    ) |> factor(levels = c("Up", "ns", "Down"))
  ) |>
  ggplot(aes(
    x     = logFC,
    y     = -log10(P.Value),
    color = direction
  )) +
  geom_point(
    size =2
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0,            linetype = "dashed", color = "grey50") +
  ggrepel::geom_label_repel(
    data = \(d) d |>
      filter(direction == "Down") |>
      slice_min(logFC, n = 5),
    aes(label = feature_id),
    color              = "#5399b0",
    fill               = "white",
    size               = 3,
    fontface           = "italic",
    box.padding        = 0.5,
    min.segment.length = 0,
    seed               = 42
  ) +
  scale_color_manual(values = c(
    "Up"   = "#9B6981FF",
    "ns"   = "grey75",
    "Down" = "#5399b0"
  )) +
  labs(
    title = "Proteomics",
    x     = expression("Log"[2]*"FC"),
    y     = expression("-Log"[10]*"P"),
    color = "Direction"
  ) +
  theme_custom()+
  theme(
    axis.text.x=element_text(angle=0,hjust=.6)
  )

### CD4 Volcano plot ------------------

cd4_plot <- da_res_total |>
  filter(grepl("CD4 T cell RNA", exp_name)) |>
  filter(!grepl("^ENSG|^LINC|-AS1$",feature_id)) |> 
  mutate(
    sig = P.Value < 0.05,
    direction = case_when(
      P.Value < 0.05 & logFC > 0 ~ "Up",
      P.Value < 0.05 & logFC < 0 ~ "Down",
      .default = "ns"
    ) |> factor(levels = c("Up", "ns", "Down"))
  ) |>
  ggplot(aes(
    x     = logFC,
    y     = -log10(P.Value),
    color = direction
  )) +
  geom_point(
    size = 2
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0,            linetype = "dashed", color = "grey50") +
  ggrepel::geom_label_repel(
    data = \(d) d |>
      filter(direction == "Down") |>
      slice_min(logFC, n = 5),
    aes(label = feature_id),
    color              = "#5399b0",
    fill               = "white",
    size               = 3,
    fontface           = "italic",
    box.padding        = 0.5,
    min.segment.length = 0,
    seed               = 42
  ) +
  scale_color_manual(values = c(
    "Up"   = "#9B6981FF",
    "ns"   = "grey75",
    "Down" = "#5399b0"
  )) +
  labs(
    title = bquote(bold("CD4")^"+" ~ bold("T cell RNA")),
    x     = expression("Log"[2]*"FC"),
    y     = expression("-Log"[10]*"P"),
    color = "Direction"
  ) +
  theme_custom()+
  theme(
    axis.text.x=element_text(angle=0,hjust=.6)
  )

### CD16 Volcano plot ------------------

cd16_plot <- da_res_total |>
  filter(grepl("CD16 Monocyte RNA", exp_name)) |>
  filter(!grepl("^ENSG|^LINC|-AS1$",feature_id)) |> 
  mutate(
    sig = P.Value < 0.05,
    direction = case_when(
      P.Value < 0.05 & logFC > 0 ~ "Up",
      P.Value < 0.05 & logFC < 0 ~ "Down",
      .default = "ns"
    ) |> factor(levels = c("Up", "ns", "Down"))
  ) |>
  ggplot(aes(
    x     = logFC,
    y     = -log10(P.Value),
    color = direction
  )) +
  geom_point(
    size =2
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0,            linetype = "dashed", color = "grey50") +
  ggrepel::geom_label_repel(
    data = \(d) d |>
      filter(direction == "Down") |>
      slice_min(logFC, n = 5),
    aes(label = feature_id),
    color              = "#5399b0",
    fill               = "white",
    size               = 3,
    fontface           = "italic",
    box.padding        = 0.5,
    min.segment.length = 0,
    seed               = 42
  ) +
  scale_color_manual(values = c(
    "Up"   = "#9B6981FF",
    "ns"   = "grey75",
    "Down" = "#5399b0"
  )) +
  labs(
    title = bquote(bold("CD16")^"+" ~ bold("Monocyte RNA")),
    x     = expression("Log"[2]*"FC"),
    y     = expression("-Log"[10]*"P"),
    color = "Direction"
  ) +
  theme_custom()+
  theme(
    axis.text.x=element_text(angle=0,hjust=.6)
  )

### Combined -------------------------
(a_plot|p_plot|cd16_plot|cd4_plot)
## Per Cytokine  Plot ------------
da_res_total |> 
  filter(exp_name == "CD4 T cell RNA") |> 
  filter(feature_id %in% c(
    msig |> filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |> pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories,
    by = c("feature_id" = "gene_symbol")
  ) |> 
  mutate(
    sig = case_when(
      P.Value < 0.05 & logFC > 0 ~ "up",
      P.Value < 0.05 & logFC < 0 ~ "down",
      .default = "ns"
    ),
    feature_id = fct_reorder(feature_id, logFC),
    category   = fct_reorder(category, logFC, .fun = mean)
  ) |> 
  ggplot(aes(
    x    = "Effect",
    y    = feature_id,
    fill = logFC
  )) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(sig != "ns", "*", "")),
    color = "white",
    size  = 3.5,
    vjust = 0.75
  ) +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey95",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Log"[2]*"FC")
  ) +
  facet_grid(
    category ~ .,
    scales = "free_y",
    space  = "free_y"
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "CD4 T cell — cytokine activity"
  ) +
  theme_custom() +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(face = "italic", size = 8),
    strip.text.y      = element_text(angle = 0, face = "bold", hjust = 0),
    panel.grid        = element_blank(),
    panel.spacing.y   = unit(0.15, "lines"),
    legend.position   = "right"
  )

### CD4 Median Cytokine ----------------
da_res |> 
  filter(exp_name == "CD4 T cell RNA") |> 
  filter(feature_id %in% c(
    msig |> filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |> pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories,
    by = c("feature_id" = "gene_symbol")
  ) |> 
  filter(!category %in% c(
    "Misc / Less Characterized",
    "Other Chemokines"
  )) |> 
  summarise(
    median_logFC = median(logFC),
    n            = dplyr::n(),
    .by          = category
  ) |> 
  mutate(
    category = fct_reorder(category, median_logFC),
    direction = ifelse(median_logFC > 0, "up", "down")
  ) |> 
  filter(n>2) |> 
  ggplot(aes(
    x    = median_logFC,
    y    = category,
    fill = median_logFC
  )) +
  geom_col(width = .85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey90",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Median Log"[2]*"FC")
  ) +
  labs(
    x     = expression("Median Log"[2]*"FC"),
    y     = NULL,
    title = bquote(bold("CD4")^"+" ~ bold("T cell"))
  ) +
  theme_custom() +
  xlim(c(-3,5.5))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    axis.text.x=element_text(angle=0,hjust=.6)
  )

### CD16 Median Cytokine ----------------
da_res |> 
  filter(exp_name == "CD16 Monocyte RNA") |> 
  filter(feature_id %in% c(
    msig |> filter(gs_name == "GOMF_CYTOKINE_ACTIVITY") |> pull(gene_symbol)
  )) |> 
  inner_join(
    cytokine_categories,
    by = c("feature_id" = "gene_symbol")
  ) |> 
  filter(!category %in% c(
    "Misc / Less Characterized",
    "Other Chemokines"
  )) |> 
  summarise(
    median_logFC = median(logFC),
    n            = dplyr::n(),
    .by          = category
  ) |> 
  mutate(
    category = fct_reorder(category, median_logFC),
    direction = ifelse(median_logFC > 0, "up", "down")
  ) |> 
  filter(n>2) |> 
  ggplot(aes(
    x    = median_logFC,
    y    = category,
    fill = median_logFC
  )) +
  geom_col(width = .85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey90",
    high     = "#9B6981FF",
    midpoint = 0,
    name     = expression("Median Log"[2]*"FC")
  ) +
  labs(
    x     = expression("Median Log"[2]*"FC"),
    y     = NULL,
    title = bquote(bold("CD16")^"+" ~ bold("Monocyte"))
  ) +
  theme_custom() +
  xlim(c(-3,5.5))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x=element_text(angle=0,hjust=.6),
    legend.position    = "none"
  )

## Save Data -----------------

## example Plots -----------------
iris |>
  filter(Petal.Length>4) |>  
  ggplot(aes(x=Sepal.Length,y=Petal.Length))+
  geom_point(color="black",
             fill="grey50",
             shape=21,
             size=1,
             alpha=.6)+
  theme_custom()+
  stat_smooth(method=lm,color="#9B6981FF")+
  # labs(x="Exposure",
  #      y="Lung Function")+
  # labs(x="Omics",
  #      y="Lung Function")+
  # labs(x="Omics",
  #      y="Exposures")+
  labs(x="Cell Types",
       y="Exposures")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
