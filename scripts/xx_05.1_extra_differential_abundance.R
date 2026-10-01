

normalize_rna_assay <- function(mae, assay_name) {
  se <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  mat <- SummarizedExperiment::assay(se, "counts") |> as.matrix()
  
  dge <- edgeR::DGEList(counts = mat) |>
    edgeR::calcNormFactors(method = "TMM")
  # 
  # keep <- edgeR::filterByExpr(dge)
  # dge  <- dge[keep, keep.lib.sizes = FALSE]
  
  log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  
  # shift to non-negative to avoid the pseudocount patch in run_differential_abundance
  # this preserves relative differences, which is all limma_trend needs
  log_cpm_shifted <- log_cpm - min(log_cpm, na.rm = TRUE)
  
  SummarizedExperiment::assay(se, "counts") <- log_cpm_shifted
  MultiAssayExperiment::experiments(mae)[[assay_name]] <- se
  mae
}

# MultiAssayExperiment Object
fev1_fvc_expom <- readRDS("./results/exwas/fev1_fvc_expom.rds")

run_deseq <- function(
    experiment_name
){
  tidyexposomics:::.update_assay_colData(
    fev1_fvc_expom,
    experiment_name)|> 
    (\(se) {
      mat <- assay(se,"counts") |>
        assay() |> 
        as.data.frame() |> 
        mutate_all(as.integer) |> 
        as.matrix()
      assay(se,"counts") <- mat
      se
    } )() |> 
    identify_abundant() |> 
    keep_abundant() |> 
    (\(se){
      res <- se |> 
        test_differential_expression(
          base_formula,
          method = "DESeq2",
        ) |> 
        (\(se) se@metadata$tidybulk$DESeq2_fit)() |> 
        mutate(exp_name = experiment_name) |> 
        inner_join(
          pivot_feature(fev1_fvc_expom),
          by=c("transcript"=".feature")
        )
    })() 
}

da_res_total <- bind_rows(
  run_deseq("CD4 T cell RNA"),
  run_deseq("CD16 Monocyte RNA"),
    tidyexposomics:::.update_assay_colData(
      fev1_fvc_expom,
      "Adductomics") |> 
      tidyexposomics:::.run_limma_trend(
        formula        = base_formula,
        abundance_col  = "counts",
        scaling_method = "none"
      ) |> 
    mutate(exp_name = "Adductomics"),
  tidyexposomics:::.update_assay_colData(
    fev1_fvc_expom,
    "Adduct Burden") |> 
    tidyexposomics:::.run_limma_trend(
      formula        = base_formula,
      abundance_col  = "counts",
      scaling_method = "none"
    ) |> 
    mutate(exp_name = "Adduct Burden"),
  tidyexposomics:::.update_assay_colData(
    fev1_fvc_expom,
    "Proteomics") |> 
    tidyexposomics:::.run_limma_trend(
      formula        = base_formula,
      abundance_col  = "counts",
      scaling_method = "none"
    ) |> 
    mutate(exp_name = "Proteomics")
) |> 
  mutate(pvalue=ifelse(
    is.na(pvalue),P.Value,pvalue
  )) |> 
  mutate(logfc=case_when(
    is.na(logFC) ~ log2FoldChange,
    is.na(log2FoldChange) ~ logFC
  ))



## No. of DEGs ----------------------
da_res_total |>
  filter(pvalue < 0.05) |> 
  dplyr::count(exp_name,sign(logfc)) |>
  mutate(direction=ifelse(
    `sign(logfc)`>0,
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
    
    "Adductomics"       = "Adductomics",
    "Adductomics Load"  = "Adductomics Load",
    "Proteomics"        = "Proteomics"
  ))+
  labs(
    x = "No. of Differential Features",
    y = "",
    fill = "Direction",
    title = bquote(bold("FEV")[1] ~ bold("/FVC")),
    subtitle = "",
    caption = "Nominal P < 0.05"
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

ranks <- da_res_total |>
  filter(!grepl("Isoform|Adductomics",exp_name)) |> 
  filter(!is.na(feature_id)) |> 
  mutate(rank_stat = sign(logfc) * -log10(pvalue)) |>
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

## Cytokine GSVA -------------------
source("./scripts/get_gsva_scores.R")
source("./scripts/pivot_se.R")
cytokine_gs_names <- msig |> 
  filter(grepl("_PRODUCTION",gs_name)) |> 
  filter(!grepl("REGULATION",gs_name)) |> 
  filter(grepl("^GOBP",gs_name)) |> 
  dplyr::count(gs_name) |> 
  pull(gs_name)


gois <- msig |> 
  filter(gs_name %in% cytokine_gs_names) |> 
  (\(df) split(df$gene_symbol,df$gs_name))()

cd4_cytokine_profiles <- get_gsva_scores(
  tidyexposomics:::.update_assay_colData(
    normalize_rna_assay(fev1_fvc_expom,"CD4 T cell RNA"),
    "CD4 T cell RNA"
  ),
  gene_sets = gois,
  gene_symbol_col = "feature_id",
  log_trans = F
)

cd16_cytokine_profiles <- get_gsva_scores(
  tidyexposomics:::.update_assay_colData(
    normalize_rna_assay(fev1_fvc_expom,"CD16 Monocyte RNA"),
    "CD16 Monocyte RNA"
  ),
  gene_sets = gois,
  gene_symbol_col = "feature_id",
  log_trans = F
)

cd4_cytokine_res <- map(names(gois),~{
  lm(
    as.formula(paste0(
      "pftfev1fvc_actual ~ ",
      .x,
      " + gli_age + gli_sex + fis"
    )),data = cd4_cytokine_profiles
  ) |> 
    broom::tidy() 
}) |> 
  bind_rows() |> 
  filter(term %in% names(gois))


cd16_cytokine_res <- map(names(gois),~{
  lm(
    as.formula(paste0(
      "pftfev1fvc_actual ~ ",
      .x,
      " + gli_age + gli_sex + fis"
    )),data = cd16_cytokine_profiles
  ) |> 
    broom::tidy() 
}) |> 
  bind_rows() |> 
  filter(term %in% names(gois))
