## Load Libraries ------------------------
library(tidyverse)
library(patchwork)
library(tidybulk)
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
cat_colors <- scales::alpha(c("#D0B8C8FF",
                              "#E08088FF", 
                              "#C8D0D8FF", 
                              "#F8C8C8FF",
                              "#A090A0FF", 
                              "#586888FF", 
                              "#182848FF",
                              "#686868FF",
                              "#F8D0B0FF", 
                              "#F8E8E8FF", 
                              "#885850FF",
                              "#B0C8F8FF",
                              "#C8E0F8FF"), 0.7) |> 
  (\(colors){
    names(colors) <- aw_cb |> 
      filter(variable %in% exp_vars) |> 
      pull(category) |> 
      unique()
    colors
  })()


exp_cols <- scales::alpha(c("#484850FF",
                            #"#000000FF", 
                            "#A0B0C0FF", 
                            #"#604830FF", 
                            "#E05060FF", 
                            "#D8B850FF",
                            "#7888A0FF",
                            "#F8F8F8FF", 
                            #"#C83850FF",
                            #"#888890FF",
                            "#883048FF",
                            #"#F8E050FF", 
                            "#D8D8D8FF",
                            "#F08070FF"
                            #"#485068FF"
                            ), 0.7) |> 
  (\(colors){
    names(colors) <- names(experiments(fev1_fvc_expom)) |> 
      unique()
    colors
  })()

msig <- msigdbr::msigdbr(species = "human",collection = "C5")

da_res_total <- readRDS("./results/da_res/da_res_total.rds")

## Function for DESeq2 --------------
run_deseq <- function(
    mae,
    experiment_name,
    fm
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
          fm,
          method = "DESeq2",
        ) |> 
        (\(se) se@metadata$tidybulk$DESeq2_fit)() |> 
        mutate(exp_name = experiment_name) |> 
        inner_join(
          pivot_feature(mae) |> 
            filter(.exp_name == experiment_name),
          by=c("transcript"=".feature")
        )
    })() 
}
## Run Exposure ~ Feature Associations -----------
exposure_formula <- function(exposure) {
  as.formula(paste("~", exposure, "+ gli_age + gli_sex + fis"))
}

da_res_exposures <- map(num_exp_vars, \(exp_var) {
  
  fm <- exposure_formula(exp_var)
  
  bind_rows(
    run_deseq(mae = fev1_fvc_expom,
              experiment_name = "CD4 T cell RNA",
              fm = fm),
    run_deseq(mae = fev1_fvc_expom,
              experiment_name = "CD16 Monocyte RNA",
              fm = fm),
    run_deseq(mae = fev1_fvc_expom,
              experiment_name = "CD4 T cell Isoforms",
              fm = fm),
    run_deseq(mae = fev1_fvc_expom,
              experiment_name = "CD16 Monocyte Isoforms",
              fm = fm),
    tidyexposomics:::.update_assay_colData(fev1_fvc_expom, "Adductomics") |>
      tidyexposomics:::.run_limma_trend(
        formula        = fm,
        abundance_col  = "counts",
        scaling_method = "none"
      ) |>
      mutate(exp_name = "Adductomics"),
    tidyexposomics:::.update_assay_colData(fev1_fvc_expom, "Adduct Burden") |>
      tidyexposomics:::.run_limma_trend(
        formula        = fm,
        abundance_col  = "counts",
        scaling_method = "none"
      ) |>
      mutate(exp_name = "Adduct Burden"),
    tidyexposomics:::.update_assay_colData(fev1_fvc_expom, "Proteomics") |>
      tidyexposomics:::.run_limma_trend(
        formula        = fm,
        abundance_col  = "counts",
        scaling_method = "none"
      ) |>
      mutate(exp_name = "Proteomics")
  ) |>
    mutate(
      exposure = exp_var,
      pvalue   = ifelse(is.na(pvalue), P.Value, pvalue),
      logfc    = case_when(
        is.na(logFC)         ~ log2FoldChange,
        is.na(log2FoldChange) ~ logFC
      )
    )
  
}) |>
  list_rbind()

## logFC Comparison ----------------
da_res_exposures <- da_res_exposures |> 
  mutate(
    unique_col =case_when(
     exp_name == "Adduct Burden"          ~ feature_id,
     exp_name == "CD16 Monocyte RNA"      ~ transcript,
     exp_name == "CD4 T cell RNA"         ~ transcript,
     exp_name == "CD16 Monocyte Isoforms" ~ transcript,
     exp_name == "CD4 T cell Isoforms"    ~ transcript,
     exp_name == "Adductomics"            ~ IonIntQuant_key,
     exp_name == "Proteomics"             ~ protein_id
    ))


da_res_total <- da_res_total |> 
  mutate(
    unique_col =case_when(
      exp_name == "Adduct Burden"          ~ feature_id,
      exp_name == "CD16 Monocyte RNA"      ~ transcript,
      exp_name == "CD4 T cell RNA"         ~ transcript,
      exp_name == "CD16 Monocyte Isoforms" ~ transcript,
      exp_name == "CD4 T cell Isoforms"    ~ transcript,
      exp_name == "Adductomics"            ~ IonIntQuant_key,
      exp_name == "Proteomics"             ~ protein_id
    ))

top_features <- da_res_exposures |>
  group_by(exp_name, unique_col) |>
  summarise(logfc_var = var(logfc, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(logfc_var), logfc_var > 0) |>
  group_by(exp_name) |>
  slice_max(logfc_var, n = 1000) |>
  ungroup()

exp_cor_res_df <- {
  names <- unique(da_res_exposures$exp_name)
  
  map(names, ~{
    da_res_exposures |>
      inner_join(top_features, by = c("exp_name", "unique_col")) |>
      filter(exp_name == .x) |>
      dplyr::select(unique_col, logfc, exposure) |>
      pivot_wider(names_from = exposure, values_from = logfc) |>
      column_to_rownames("unique_col") |>
      as.matrix() |>
      (\(m) m[complete.cases(m), ])() |>
      Hmisc::rcorr() |>
      (\(lst){
        lst$r |>
          as.data.frame() |>
          rownames_to_column("var1") |>
          pivot_longer(!var1, names_to = "var2", values_to = "r") |>
          inner_join(
            lst$P |>
              as.data.frame() |>
              rownames_to_column("var1") |>
              pivot_longer(!var1, names_to = "var2", values_to = "pvalue"),
            by = c("var1", "var2")
          )
      })() |>
      filter(var1 != var2) |>
      mutate(exp_name = .x)
  }) |>
    bind_rows()
} |>
  left_join(
    aw_cb |>
      dplyr::select(var1 = variable, var1_clean = clean_name, var1_cat = category),
    by = "var1"
  ) |>
  left_join(
    aw_cb |>
      dplyr::select(var2 = variable, var2_clean = clean_name, var2_cat = category),
    by = "var2"
  ) |>
  mutate(cat_comparison = case_when(
    var1_cat == var2_cat ~ "Within Category Comparison",
    .default = "Across Category Comparison"
  ))

### -----------

exp_cor_sum_df <- exp_cor_res_df |> 
  mutate(
    combo = map2_chr(var1_cat, 
                     var2_cat,
                     \(a, b) paste(sort(c(a, b)),
                                   collapse = "_"))
  ) |> 
  group_by(exp_name, combo) |>
  reframe(
    med_cor = median(r),
    med_logp = median(-log10(pvalue + 1e-300))
  ) |> 
  separate(combo,sep="_",
           into = c("var1_cat","var2_cat"),
           remove = F) |> 
  mutate(cat_comparison=case_when(
    var1_cat == var2_cat ~ "Within Category Comparison",
    .default = "Across Category Comparison"
  )) 
  

exp_cor_sum_df |> 
  group_by(exp_name) |> 
  arrange(desc(med_cor)) |> 
  ggplot(aes(
    x=var1_cat,
    y=var2_cat,
    fill = med_cor
  ))+
  geom_tile()+
  theme_custom()+
  theme(
    axis.text.x=element_text(angle=90,vjust=.6)
  )+
  facet_wrap(~exp_name)+
  scale_fill_gradient2(
    high  = "#C33764",
    mid   = "white",
    low   = "#1D2671",
    midpoint = 0,
    limits=c(-1,1)
  )+
  labs(
    x="",
    y="",
    fill="Median Correlation"
  )

## GSEA -----------
source("./scripts/gsea_go_enrichment.R")
exp_ranks <- da_res_exposures |>
  filter(!grepl("Isoform|Adductomics",exp_name)) |> 
  filter(!is.na(feature_id)) |> 
  mutate(rank_stat = sign(logfc) * -log10(pvalue)) |>
  filter(
    is.finite(rank_stat),
    !is.na(rank_stat)
  ) |>
  mutate(combo=paste(exp_name,exposure,sep="-")) |> 
  (\(df) split(df, df$combo))() |>
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

exp_gsea_res_sig <- map(
  exp_ranks,
  \(ranks) run_fgsea(
    ranked_genes = ranks,
    db           = "GO",
    term_data    = go_term_data,
    species      = "goa_human",
    feature_col  = "gene_symbol"
  ) |>
    filter(size>5,
           padj<0.1)
) |>
  dplyr::bind_rows(.id = "group")  |> 
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", ")) 

# saveRDS(exp_gsea_res_sig ,file = "./results/exposure_omics/exp_gsea_res_sig.rds")

exp_gsea_res_sig_clean <- exp_gsea_res_sig |> 
  separate(group,
           into = c("exp_name","exposure"),
           sep = "-",
           remove = F) |> 
  left_join(aw_cb,
            by=c("exposure"="variable"))
## No. of Omics Associations per exposure -----------------
exp_omics_assoc_df <- da_res_exposures |> 
  filter(pvalue<0.05) |> 
  dplyr::count(exposure,exp_name) |>
  inner_join(
    fev1_fvc_expom |> 
      extract_results("codebook"),
    by=c("exposure"="variable")
  ) |> 
  group_by(clean_name) |> 
  mutate(total=sum(n)) |> 
  ungroup() |> 
  # Order categories by median effect size
  mutate(
    clean_name = forcats::fct_reorder(clean_name, total),
    category = forcats::fct_reorder(
      category, total, .fun = median, .desc = TRUE),
  )
  
exp_omics_assoc_df |> 
  ggplot(aes(x = n, 
             y = clean_name,
             fill = exp_name)) +
  geom_col()+
  ggh4x::facet_grid2(category ~ ., 
                     scales = "free",
                     space = "free", 
                     strip = ggh4x::strip_themed(
                       background_y = ggh4x::elem_list_rect(
                         fill = cat_colors[levels(exp_omics_assoc_df$category)] |> 
                           scales::alpha(.3)
                       ))) +
  scale_fill_manual(values = exp_cols)+
  labs(
    title = NULL,
    x = "No. of Associations",
    y = "",
    fill ="Omics Layer"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        #legend.direction = "vertical",
        legend.title.position = "top",
        strip.text.y = element_text(angle=0,
                                    face="bold.italic"))

## Exposure Category by Omics Layer Rank Heatmap -------
# no. of deges per exposure, then get the median number of degs
# per omics layer and category
# then normalize by the number of variables in the exposure category
# then group by the omics layer and rank which exposure categories
# impact that layer the most
da_res_exposures |>
  filter(pvalue < 0.05) |>
  dplyr::count(exposure, exp_name) |>
  inner_join(
    fev1_fvc_expom |> extract_results("codebook"),
    by = c("exposure" = "variable")
  ) |>
  group_by(exp_name, category) |>
  reframe(
    med_n        = median(n),
    n_exposures  = n_distinct(exposure)
  ) |>
  mutate(norm_n = med_n / n_exposures) |>
  arrange(desc(norm_n)) |>
  group_by(exp_name) |>
  mutate(rank = percent_rank(norm_n)) |> 
  ggplot(aes(
    x=reorder(exp_name,rank),
    y=reorder(category,rank),
    fill=rank
  ))+
  geom_tile()+
  theme_custom()+
  scale_fill_gradientn(
    colors = c("#0B0405FF","#357BA2FF","#DEF5E5FF") |> scales::alpha(.8) |>  rev()
  )+
  scale_x_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
    "CD4 T cell Isoforms"       = expression(CD4^"+" ~ "T cell Isoforms"),
    "CD16 Monocyte Isoforms"  = expression(CD16^"+" ~ "Monocyte Isoforms"),
    "Proteomics" ~ "Proteomics",
    "Adductomics" ~ "Adductomics",
    "Adduct Burden" ~ "Adduct Burden"
  ))+
  labs(
    x="",
    y="",
    fill="Median Rank"
  )

## Similarity with Lung Features ----------
source("./scripts/get_pairwise_overlaps.R")

jc_res_lung_df <- da_res_exposures |> 
  filter(pvalue<0.05) |> 
  mutate(combo=paste(exp_name,unique_col,sep="-")) |> 
  (\(df) split(df$combo,df$exposure))() |> 
  c(
    list("Lung Function" = da_res_total |> 
           mutate(
             unique_col =case_when(
               exp_name == "Adduct Burden"          ~ feature_id,
               exp_name == "CD16 Monocyte RNA"      ~ transcript,
               exp_name == "CD4 T cell RNA"         ~ transcript,
               exp_name == "CD16 Monocyte Isoforms" ~ transcript,
               exp_name == "CD4 T cell Isoforms"    ~ transcript,
               exp_name == "Adductomics"            ~ IonIntQuant_key,
               exp_name == "Proteomics"             ~ protein_id
             )) |> 
           filter(pvalue<0.05) |> 
           mutate(combo=paste(exp_name,unique_col,sep="-")) |> 
           pull(combo))
  ) |> 
  get_pairwise_overlaps() |> 
  filter(source == "Lung Function" |
           target == "Lung Function") |> 
  left_join(
    aw_cb,
    by=c("source"="variable")
  )

jc_res_lung_df |> 
  ggplot(aes(
    x=reorder(category,jaccard),
    y=jaccard,
    fill=category,
    color=category
  ))+
  geom_jitter(alpha=.1)+
  geom_boxplot(color="black")+
  theme_custom()+
  scale_fill_manual(values=cat_colors)+
  scale_color_manual(values = cat_colors)+
  theme(legend.position = "none")+
  labs(
    x="",
    y="Jaccard Coef."
  )
  
## Lung Function Direction Concordance -------------
lung_dir_concord <- map(unique(da_res_exposures$exposure),~{
  
  da_res_exposures |> 
    filter(exposure==.x) |> 
    mutate(direction=ifelse(logfc>0,"up","down")) |> 
    dplyr::select(unique_col,
                  exp_name,
                  exposure_direction=direction) |> 
    inner_join(
      da_res_total |> 
        mutate(direction=ifelse(logfc>0,"up","down")) |> 
        dplyr::select(unique_col,
                      exp_name,
                      lung_function_direction=direction),
      by=c("unique_col",
           "exp_name")
    ) |> 
    mutate(concordance=case_when(
      exposure_direction != lung_function_direction ~ "concordant",
      .default = "not concordant"
    )) |> 
    dplyr::count(concordance) |> 
    mutate(total=sum(n)) |> 
    filter(concordance != "not concordant") |> 
    mutate(pct_concordant=n/total) |> 
    mutate(exposure=.x)
}) |> 
  bind_rows() |> 
  left_join(
    aw_cb,
    by=c("exposure"="variable")
  )

lung_dir_concord |> 
  ggplot(aes(
    x=reorder(category,pct_concordant),
    y=pct_concordant,
    fill=category,
    color=category
  ))+
  geom_jitter(alpha=.1)+
  geom_boxplot(color="black")+
  theme_custom()+
  scale_fill_manual(values=cat_colors)+
  scale_color_manual(values = cat_colors)+
  theme(legend.position = "none")+
  scale_y_continuous(labels = scales::percent)+
  labs(
    x="",
    y="% Concordance"
  )

## DEG Concordance v. Jaccard -------------
lung_dir_concord |> 
  dplyr::select(pct_concordant,
                exposure,category) |> 
  left_join(jc_res_lung_df |> 
              dplyr::select(exposure=source,category,jaccard),
            by=c("exposure","category")) |> 
  mutate(pct_concordant_rank = percent_rank(pct_concordant),
         jaccard_rank = percent_rank(jaccard)) |> 
  group_by(category) |> 
  reframe(
    med_concord_rank = median(pct_concordant_rank),
    med_jaccard_rank = median(jaccard_rank)
  ) |> 
  ggplot(aes(
    x=med_concord_rank,
    y=med_jaccard_rank,
    fill = category
  ))+
  geom_hline(yintercept = .5,linetype="dashed")+
  geom_vline(xintercept = .5,linetype="dashed")+
  geom_point(size=5,shape=21)+
  scale_fill_manual(values=cat_colors)+
  theme_custom()+
  labs(
    x="Median % Concordance Rank",
    y="Median Jaccard Coef. Rank",
    fill="Category",
    title="Omics Features"
  )

## Pathway Concordance v. Jaccard --------
exp_gsea_res_sig_clean |> 
  inner_join(
    gsea_res_sig |> 
      dplyr::select(
        group,
        term_name,
        lung_nes = NES),
    by=c("exp_name"="group",
         "term_name"="term_name")
  ) |> 
  mutate(concordant=sign(NES) != sign(lung_nes)) |>
  group_by(exposure) |> 
  reframe(n_concord=length(concordant[concordant==TRUE])) |> 
  inner_join(
    exp_gsea_res_sig_clean |>
      dplyr::count(exposure) |> 
      dplyr::select(
        exposure,
        total_path_n=n
      ),
    by="exposure"
  ) |> 
  mutate(pct_concord=n_concord/total_path_n) |> 
  inner_join(
    exp_gsea_res_sig_clean |> 
      (\(df) split(df$term_name,df$exposure))() |> 
      c(list(
        "Lung Function"=gsea_res_sig |> pull(term_name)
      )) |> 
      map(~{
        .x |> unique()
      }) |> 
      get_pairwise_overlaps() |> 
      filter(source == "Lung Function" |
               target == "Lung Function") |> 
      left_join(
        aw_cb,
        by=c("source"="variable")
      ),
    by=c("exposure"="source")
  ) |> 
  mutate(
    rank_concord=percent_rank(pct_concord),
    rank_jaccard=percent_rank(jaccard)
  ) |> 
  group_by(category) |> 
  reframe(
    med_concord_rank = median(rank_concord),
    med_jaccard_rank = median(rank_jaccard)
  ) |> 
  ggplot(aes(
    x=med_concord_rank,
    y=med_jaccard_rank,
    fill = category
  ))+
  geom_hline(yintercept = .5,linetype="dashed")+
  geom_vline(xintercept = .5,linetype="dashed")+
  geom_point(size=5,shape=21)+
  scale_fill_manual(values=cat_colors)+
  theme_custom()+
  labs(
    x="Median % Concordance Rank",
    y="Median Jaccard Coef. Rank",
    fill="Category",
    title="Pathways"
  )

## Which omics layers are concordant ------------
exp_gsea_res_sig_clean |> filter(size<200) |> 
  inner_join(
    gsea_res_sig |>
      dplyr::select(group, term_name, lung_nes = NES),
    by = c("exp_name" = "group", "term_name" = "term_name")
  ) |>
  mutate(concordant = sign(NES) != sign(lung_nes)) |>
  filter(concordant == TRUE) |> 
  dplyr::count(exp_name) 

## Which pathways are being affected? --------
exp_gsea_res_sig_clean |> 
  inner_join(
    gsea_res_sig |>
      dplyr::select(group, term_name, lung_nes = NES),
    by = c("exp_name" = "group", "term_name" = "term_name")
  ) |>
  mutate(concordant = sign(NES) != sign(lung_nes)) |>
  filter(concordant == TRUE) |> 
  dplyr::count(term_name) |>
  arrange(desc(n)) |> 
  # removing super broad terms
  filter(!term_name %in% c(
    "cytosol","nucleus","cytoplasm"
  )) |> 
  slice_head(n=10) |> 
  ggplot(aes(
    x=n,
    y=reorder(term_name,n),
    fill=n
  ))+
  geom_col(width=.9)+
  theme_custom()+
  scale_fill_gradientn(
    colors = c("#000004FF","#B63679FF","#FCFDBFFF") |> 
      rev() |> 
      scales::alpha(.8),
    
  )+
  theme(
    axis.text.x = element_text(angle=0,hjust = 0.2)
  )+
  labs(
    x="No. of Exposures",
    y="",
    fill="No. of Exposures"
  )
  
## Which exposure categories are driving this --------
exp_gsea_res_sig_clean |> 
  inner_join(
    gsea_res_sig |>
      dplyr::select(group, term_name, lung_nes = NES),
    by = c("exp_name" = "group", "term_name" = "term_name")
  ) |>
  mutate(concordant = sign(NES) != sign(lung_nes)) |>
  filter(concordant == TRUE) |> 
  # removing super broad terms
  filter(!term_name %in% c(
    "cytosol","nucleus","cytoplasm"
  )) |> 
  group_by(category) |> 
  reframe(n_path=n_distinct(term_name)) |> 
  arrange(desc(n_path)) |>  
  slice_head(n=10) |> 
  ggplot(aes(
    x=n_path,
    y=reorder(category,n_path),
    fill=category
  ))+
  geom_col(width=.85)+
  theme_custom()+
  scale_fill_manual(values = cat_colors)+
  theme(
    axis.text.x = element_text(angle=0,hjust = 0.2)
  )+
  guides(fill="none")+
  labs(
    x="No. of Pathways",
    y=""
  )



## Top Lung Concordant Pathways across exposures ----------


exp_gsea_res_sig_clean |>
  inner_join(
    gsea_res_sig |>
      dplyr::select(group, term_name, lung_nes = NES),
    by = c("exp_name" = "group", "term_name" = "term_name")
  ) |>
  mutate(concordant = sign(NES) != sign(lung_nes)) |>
  filter(concordant == TRUE) |>
  (\(d) {
    cat_counts <- d |>
      group_by(term_name, exp_name, category) |>
      reframe(n_cat_exposures = n_distinct(exposure))
    
    d |>
      group_by(term_name, exp_name) |>
      reframe(
        n_exposures = n_distinct(exposure),
        categories  = paste(unique(category), collapse = ","),
        med_nes     = median(NES)
      ) |>
      mutate(
        top5 = term_name %in% (pick(everything()) |>
                                 group_by(exp_name) |>
                                 slice_max(n_exposures,
                                           n = 7,
                                           with_ties = FALSE) |>
                                 pull(term_name))
      ) |>
      filter(top5) |>
      mutate(term_name = fct_reorder(term_name, n_exposures, mean)) |>
      (\(d2) {
        term_levels <- levels(d2$term_name)
        
        p_nes <- d2 |>
          mutate(term_name = factor(term_name, levels = term_levels)) |>
          ggplot(aes(
            x = exp_name,
            y = term_name,
            fill = med_nes)) +
          geom_tile() +
          scale_fill_gradient2(
            low      = "#5399b0",
            mid      = "grey85",
            high     = "#9B6981FF",
            midpoint = 0,
            name     = "Median NES"
          ) +
          labs(x = NULL, 
               y = NULL, 
               title = "Median NES") +
          theme_custom() +
          theme(axis.text.x = element_text(angle = 35, hjust = 1),
                legend.title.position = "top")
        
        p_bar <- cat_counts |>
          filter(term_name %in% term_levels) |>
          mutate(term_name = factor(term_name,
                                    levels = term_levels)) |>
          group_by(term_name, category) |>
          reframe(n_cat_exposures = sum(n_cat_exposures)) |>
          ggplot(aes(x = n_cat_exposures, y = term_name, fill = category)) +
          geom_col(position = "fill") +
          scale_x_continuous(labels = scales::percent) +
          scale_fill_manual(values = cat_colors, 
                            name = "Exposure category") +
          labs(x = "% Exposures",
               y = NULL, 
               title = "Category breakdown") +
          theme_custom() +
          theme(axis.text.y = element_blank())
        
        p_dot <- d2 |>
          mutate(term_name = factor(term_name,
                                    levels = term_levels)) |>
          group_by(term_name) |>
          reframe(n_exposures = sum(n_exposures)) |>
          ggplot(aes(x = n_exposures, y = term_name)) +
          geom_segment(
            aes(x = 0, xend = n_exposures, yend = term_name),
            colour = "grey70",
            linewidth = 0.4
          ) +
          geom_point(size = 3, shape = 21, fill = "grey30") +
          labs(x = "No. exposures", 
               y = NULL, 
               title = "Total exposures") +
          theme_custom() +
          theme(axis.text.y = element_blank())
        
        p_nes + p_bar + p_dot +
          plot_layout(widths = c(1, 3, 2), guides = "collect") &
          theme(legend.position = "bottom")
      })()
  })()
## Similarity between exposure-omics features --------------
source("./scripts/get_pairwise_overlaps.R")

jc_res_df <- da_res_exposures |> 
  filter(pvalue<0.05) |> 
  mutate(combo=paste(exp_name,unique_col,sep="-")) |> 
  (\(df) split(df$combo,df$exposure))() |> 
  get_pairwise_overlaps() |> 
  left_join(
    aw_cb |>
      dplyr::select(source = variable, 
                    source_clean = clean_name,
                    source_cat = category),
    by = "source"
  ) |>
  left_join(
    aw_cb |>
      dplyr::select(target = variable, 
                    target_clean = clean_name, 
                    target_cat = category),
    by = "target"
  ) |> 
  mutate(cat_comparison = case_when(
    source_cat == target_cat ~ "Within Category Comparison",
    .default = "Across Category Comparison"
  )) |> 
  view()

### JC Per Exposure Plotting ---------------
jc_sum_df <- jc_res_df |> 
  dplyr::select(source,target,jaccard) |> 
  pivot_longer(
    cols = source:target,
    names_to = "exposure_cat",
    values_to = "exposure"
  ) |> 
  group_by(exposure) |> 
  reframe(
    med_jc=median(jaccard)
  ) |> 
  inner_join(
    aw_cb,
    by=c("exposure"="variable")
  ) |> 
  mutate(
    clean_name = tidytext::reorder_within(
      x = clean_name,
      by = med_jc,
      within = category
    ),
    category = forcats::fct_reorder(
      category, med_jc, .fun = median, .desc = TRUE)) 


jc_sum_df |>
  ggplot(aes(x = med_jc, 
             y = clean_name,
             color = med_jc,
             fill = med_jc)) +
  geom_point(shape = 21, 
             size = 2.5, 
             stroke = 0.8) +
  geom_segment(aes(
    x = 0,
    xend = med_jc,
    y=clean_name
  ))+
  tidytext::scale_y_reordered()+
  ggh4x::facet_grid2(category ~ ., 
                     scales = "free",
                     space = "free", 
                     strip = ggh4x::strip_themed(
                       background_y = ggh4x::elem_list_rect(
                         fill = cat_colors[levels(jc_sum_df$category)] |> 
                           scales::alpha(.3)
                         
                       ))) +
  scale_color_gradientn(
    colors = c("#000004FF","#B63679FF","#FCFDBFFF") |> rev()
  )+
  scale_fill_gradientn(
    colors = c("#000004FF","#B63679FF","#FCFDBFFF") |> rev(),
    
  )+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = NULL,
    x = "Median JC",
    y = "",
    fill ="Median JC",
    color ="Median JC"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title.position = "top",
        strip.text.y = element_text(angle=0,
                                    face="bold.italic"))

### JC v. DEGs ----------
jc_sum_df |> 
  mutate(clean_name=gsub("__.*","",clean_name)) |> 
  inner_join(
    da_res_exposures |>
      filter(pvalue < 0.05) |>
      dplyr::count(exposure, name = "n_deg"),
    by = "exposure"
  ) |>
  mutate(
    label = if_else((n_deg > 3000 & med_jc < 0.5), clean_name, NA)
  ) |>
  ggplot(aes(
    x     = n_deg,
    y     = med_jc,
    color = category,
    label = label
  )) +
  geom_point(size = 3) +
  ggrepel::geom_label_repel(size = 3.5,
                            max.overlaps = 30,
                            force = 30,
                            color = "black", 
                            fill="white",
                            na.rm = TRUE) +
  scale_color_manual(values = cat_colors) +
  labs(
    x     = expression("No. of DEGs"),
    y     = "Median JC",
    color = "Category"
  ) +
  guides(color = guide_legend(nrow = 6)) +
  theme_custom()+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top")
## Exposure -omics Association Rank ------------

# which exposures are consistently towards the top per omics layer as having the highest number of differentially expressed features?
da_res_exposures |>
  filter(pvalue < 0.05) |>
  dplyr::count(exposure, exp_name) |>
  inner_join(
    fev1_fvc_expom |> extract_results("codebook"),
    by = c("exposure" = "variable")
  ) |>
  group_by(exp_name, category) |>
  reframe(
    med_n        = median(n),
    n_exposures  = n_distinct(exposure)
  ) |>
  mutate(norm_n = med_n / n_exposures) |>
  arrange(desc(norm_n)) |>
  group_by(exp_name) |>
  mutate(rank = percent_rank(norm_n)) |>
  group_by(category) |>
  reframe(med_rank = median(rank)) |>
  ggplot(aes(
    x    = med_rank,
    y    = reorder(category, med_rank),
    fill = category
  )) +
  geom_col() +
  theme_custom() +
  scale_fill_manual(values = cat_colors) +
  theme(legend.position = "none") +
  labs(
    x = "Median Rank",
    y = ""
  )

## No. of Common DEGs -------------------
da_res_exposures |> 
  filter(pvalue < 0.05) |> 
  dplyr::select(exposure,
                unique_col,
                exp_name) |> 
  mutate(combo=paste(unique_col,exp_name,sep=":")) |> 
  left_join(aw_cb,
            by=c("exposure"="variable")) |> 
  dplyr::select(combo,category) |> 
  distinct() |> 
  dplyr::count(combo) |> 
  filter(n>6) |> 
  separate(combo,
           into = c("unique_col","exp_name"),
           remove = F,
           sep=":") |> 
  inner_join(
    da_res_exposures |> 
      dplyr::select(unique_col,
                    exp_name,
                    feature_id,
                    feature_map),
    by=c("exp_name","unique_col")
  ) |> 
  distinct() |> 
  dplyr::count(exp_name) |> 
  ggplot(aes(
    x=n,
    y=reorder(exp_name,n),
    fill = n
  ))+
  geom_col()+
  theme_custom()+
  scale_fill_gradientn(
    colors = c("#000004FF","#B63679FF","#FCFDBFFF") |> rev(),
    
  )+
  labs(
    x="No. of DEGs",
    y="",
    fill="No. of DEGs"
  )
## GO for common DEGs -------------
source("./scripts/gsea_go_enrichment.R")
go_term_data <- read_rds("./results/useful_data/go_term_data.rds")

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

exp_omic_da_lst <- da_res_exposures |> 
  filter(pvalue < 0.05) |> 
  dplyr::select(exposure,
                unique_col,
                exp_name) |> 
  mutate(combo=paste(unique_col,exp_name,sep=":")) |> 
  left_join(aw_cb,
            by=c("exposure"="variable")) |> 
  dplyr::select(combo,category) |> 
  distinct() |> 
  dplyr::count(combo) |> 
  filter(n>6) |> 
  separate(combo,
           into = c("unique_col","exp_name"),
           remove = F,
           sep=":") |> 
  inner_join(
    da_res_exposures |> 
      dplyr::select(unique_col,
                    exp_name,
                    feature_id,
                    feature_map),
    by=c("exp_name","unique_col")
  ) |> 
  distinct() |> 
  (\(df) split(df$feature_id,df$exp_name))() |> 
  (\(lst) lst[lengths(lst) > 10])()

# Build universes per cell type by modality
universes <- pivot_feature(fev1_fvc_expom) |>
  filter(.exp_name %in% c(
    "CD4 T cell RNA", "CD4 T cell Isoforms",
    "CD16 Monocyte RNA", "CD16 Monocyte Isoforms"
  )) |>
  dplyr::count(.exp_name, feature_id) |>   # deduplicates
  split(~.exp_name) |>
  map(~ pull(.x, feature_id))

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
exp_omic_enr_res <- future_imap(
  exp_omic_da_lst,
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

exp_omic_enr_res_sig <- exp_omic_enr_res |>
  filter(
    p_adjust < .1,
    n_with_sel > 5,
    N_with < 1000
  )

### GO Plotting ------------------
exp_omic_enr_res_sig |>
  mutate(group=gsub("_","-",group)) |> 
  group_by(group) |>
  arrange(desc(-log10(p_adjust))) |>
  slice_head(n = 15) |>
  ungroup() |>
  mutate(term_label = str_trunc(term_name, width = 40)) |>
  mutate(term_label = reorder(
    term_label,
    as.integer(factor(group)) * 1e6 + -log10(p_adjust))) |>
  ggplot(aes(
    x = group,
    y = term_label,
    fill = -log10(p_adjust)
  )) +
  geom_tile() +
  theme_custom() +
  scale_fill_gradient(
    high = "#360033",
    low  = "#FFFDE4"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .4)
  ) +
  labs(
    x = "",
    y = "",
    fill = expression("-Log"[10]*"P")
  )
### Focused Sets ------------------

#### Interferon --------------
da_res_exposures |> 
  filter(pvalue<0.05) |> 
  filter(grepl("RNA",exp_name)) |> 
  filter(feature_id %in% c(
    exp_omic_enr_res_sig |> 
      filter(grepl("interferon",term_name)) |> 
      separate_rows(ids,sep=", ") |> 
      pull(ids)
  )) |> 
  group_by(exp_name,feature_id) |> 
  reframe(
    n_exposures=n_distinct(exposure),
    n_directions=n_distinct(sign(logfc)),
    med_logfc=median(logfc),
    med_logp=median(-log10(pvalue))
  ) |> 
  arrange(desc(abs(med_logfc)),desc(n_exposures)) |> 
  view()


da_res_exposures |> 
  filter(pvalue<0.05) |> 
  filter(grepl("RNA",exp_name)) |> 
  filter(feature_id %in% c(
    exp_omic_enr_res_sig |> 
      filter(grepl("chromatin|histone",term_name)) |> 
      separate_rows(ids,sep=", ") |> 
      pull(ids)
  )) |> 
  group_by(exp_name) |> 
  reframe(
    n_exposures=n_distinct(exposure),
    n_directions=n_distinct(sign(logfc)),
    med_logfc=median(logfc),
    med_logp=median(-log10(pvalue))
  ) |> 
  view()

## Discordance Analysis -------------
discordant_exposures <- jc_sum_df |>
  mutate(clean_name = gsub("__.*", "", clean_name)) |>
  inner_join(
    da_res_exposures |>
      filter(pvalue < 0.05) |>
      dplyr::count(exposure, name = "n_deg"),
    by = "exposure"
  ) |>
  filter(n_deg > 3000 & med_jc < 0.5) |>
  dplyr::pull(exposure)

# Features that are frequently DA across exposures (the "common" set to exclude)
frequent_combos <- da_res_exposures |>
  filter(pvalue < 0.05) |>
  dplyr::select(exposure, unique_col, exp_name) |>
  mutate(combo = paste(unique_col, exp_name, sep = ":")) |>
  left_join(aw_cb, by = c("exposure" = "variable")) |>
  dplyr::select(combo, category) |>
  distinct() |>
  dplyr::count(combo) |>
  filter(n > 6) |>
  dplyr::pull(combo)

# Discordant exposure DEGs, excluding frequent set,
# keeping only features seen in >= 2 discordant exposures
discordant_da_lst <- da_res_exposures |>
  filter(pvalue < 0.05, exposure %in% discordant_exposures) |>
  dplyr::select(exposure, unique_col, exp_name, feature_id) |>
  mutate(combo = paste(unique_col, exp_name, sep = ":")) |>
  filter(!combo %in% frequent_combos) |>
  distinct() |>
  add_count(unique_col, exp_name, name = "n_discordant_exposures") |>
  filter(n_discordant_exposures >= 2) |>
  mutate(group = paste(exposure, exp_name, sep = " | ")) |>
  (\(df) split(df$feature_id, df$group))() |>
  (\(lst) lst[lengths(lst) > 10])()

discordant_enr_res <- future_imap(
  discordant_da_lst,
  ~ {
    options(GO_ANNOTATION_URL = "http://current.geneontology.org/annotations/goa_human.gaf.gz")
    run_fenr(
      selected_genes = .x,
      term_data      = go_term_data,
      universe_genes = pick_universe(.y),
      db             = "GO",
      species        = "goa_human",
      feature_col    = "gene_symbol"
    )
  },
  .progress = TRUE
) |>
  bind_rows(.id = "group")

discordant_enr_res_sig <- discordant_enr_res |>
  filter(
    p_adjust   < 0.1,
    n_with_sel > 5,
    N_with     < 1000
  ) |> 
  separate(group,sep="\\| ",into = c("exposure","exp_name"),remove = F) |> 
  mutate(exposure=str_trim(exposure)) |> 
  mutate(exp_name=str_trim(exp_name))





discordant_enr_res_sig |>
  group_by(group) |>
  arrange(desc(-log10(p_adjust))) |>
  slice_head(n = 5) |>
  ungroup() |>
  inner_join(aw_cb,
             by=c("exposure" ="variable")) |> 
  mutate(term_label = str_trunc(term_name, width = 40)) |>
  mutate(term_label = reorder(
    term_label,
    as.integer(factor(group)) * 1e6 + -log10(p_adjust))) |>
  mutate(term_label = tidytext::reorder_within(term_label,by = -log10(p_adjust),within = clean_name)) |> 
  ggplot(aes(
    x = exp_name,
    y = term_label,
    fill = -log10(p_adjust)
  )) +
  geom_tile() +
  tidytext::scale_y_reordered()+
  theme_custom() +
  scale_fill_gradient(
    high = "#360033",
    low  = "#FFFDE4"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .4)
  ) +
  facet_wrap(~clean_name,nrow=4,space="free_y",scales="free_y")+
  labs(
    x = "",
    y = "",
    fill = expression("-Log"[10]*"P")
  )

### Overlap between discordant and common terms ---------

list("Discordant Terms"=discordant_enr_res_sig$term_name |> unique(),
     "Concordant Terms"=exp_omic_enr_res_sig$term_name |> unique()) |> 
  ggvenn::ggvenn(
    stroke_size = .5,
    fill_color = c("#9D1B1FFF","#A8CDECFF"),set_name_size = 5)

## Exposure-Pathway Rank ---------------
exp_gsea_res_sig_clean |>
  group_by(exp_name, category) |>
  reframe(
    n_terms     = n_distinct(term_name),
    med_nes     = median(NES),
    mean_abs_nes = mean(abs(NES)),
    pct_up      = mean(NES > 0)
  ) |>
  inner_join(
    aw_cb |> filter(variable %in% num_exp_vars) |> dplyr::count(category),
    by = "category"
  ) |>
  mutate(norm_path_count = n_terms / n) |>
  group_by(exp_name) |>
  mutate(rank = percent_rank(norm_path_count)) |>
  ungroup() |>
  group_by(category) |>
  reframe(
    med_rank      = median(rank),
    med_nes       = median(med_nes),       # median of per-omic medians
    med_abs_nes   = median(mean_abs_nes),  # effect size regardless of direction
    med_pct_up    = median(pct_up)         # >0.5 = mostly upregulated across omics
  ) |>
  arrange(desc(med_rank)) |> 
  ggplot(aes(
    x    = med_rank,
    y    = reorder(category, med_rank),
    fill = category
  )) +
  geom_col() +
  theme_custom() +
  scale_fill_manual(values = cat_colors) +
  theme(legend.position = "none") +
  labs(
    x = "Median Rank",
    y = ""
  )

## Isoform Impact --------------



iso_only_res <- map(unique(da_res_exposures$exposure),~{
  
  cd4_res <- da_res_exposures |> 
    filter(pvalue<0.05) |> 
    filter(exposure == .x) |> 
    filter(grepl("CD4",exp_name)) |> 
    (\(df) split(df$feature_id,df$exp_name))() |> 
    map(~{.x |> unique()})
  
  cd4_iso_only <- cd4_res$`CD4 T cell Isoforms`[
    !cd4_res$`CD4 T cell Isoforms` %in% cd4_res$`CD4 T cell RNA`] |> 
    length()
  
  cd4_total_hits <- n_distinct(da_res_exposures |> 
                                 filter(pvalue<0.05) |> 
                                 filter(exposure == .x) |> 
                                 filter(grepl("CD4",exp_name)) |>
                                 pull(feature_id))
  
  
  cd16_res <- da_res_exposures |> 
    filter(pvalue<0.05) |> 
    filter(exposure == .x) |> 
    filter(grepl("CD16",exp_name)) |> 
    (\(df) split(df$feature_id,df$exp_name))() |> 
    map(~{.x |> unique()})
  
  cd16_iso_only <- cd16_res$`CD16 Monocyte Isoforms`[
    !cd16_res$`CD16 Monocyte Isoforms` %in% cd16_res$`CD16 Monocyte RNA`] |> 
    length()
  
  cd16_total_hits <- n_distinct(da_res_exposures |> 
                                  filter(pvalue<0.05) |> 
                                  filter(exposure == .x) |> 
                                  filter(grepl("CD16",exp_name)) |>
                                  pull(feature_id))
  
  res <- data.frame(
    cd4_iso_only = cd4_iso_only,
    cd4_total_hits = cd4_total_hits,
    cd4_iso_only_pct = cd4_iso_only/cd4_total_hits,
    cd16_iso_only = cd16_iso_only,
    cd16_total_hits = cd16_total_hits,
    cd16_iso_only_pct = cd16_iso_only/cd16_total_hits,
    exposure=.x
  )
  
}) |> 
  bind_rows() |> 
  left_join(
    aw_cb |> 
      filter(variable %in% num_exp_vars),
    by=c("exposure"="variable")
  )


iso_only_res |> 
  group_by(category) |> 
  reframe(
    med_cd4_pct = median(cd4_iso_only_pct),
    med_cd16_pct = median(cd16_iso_only_pct)
  ) |> 
  ggplot(aes(
    x=med_cd16_pct,
    y=med_cd4_pct,
    fill=category
  ))+
  geom_point(shape=21,size=6)+
  scale_fill_manual(values=cat_colors)+
  theme_custom()+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  labs(
    x=expression("Median CD16"^"+" ~ "Monocyte Isoform Only %"),
    y=expression("Median CD4"^"+" ~ "T cell Isoform Only %"),
    fill="Category"
  )
## Old Isoform Impact --------------
isoform_assays <- c("CD16 Monocyte Isoforms", "CD4 T Cell Isoforms")
rna_assays     <- c("CD16 Monocyte RNA", "CD4 T Cell RNA")  # paired gene-level assays

# Significant hits at gene-level RNA (to exclude)
gene_sig <- da_res_exposures |>
  filter(exp_name %in% rna_assays, pvalue < 0.05) |>
  dplyr::select(exposure, feature_id) |>
  distinct()

# Isoform-only: significant in isoform assay but NOT significant at gene-level RNA
isoform_only <- da_res_exposures |>
  filter(exp_name %in% isoform_assays, pvalue < 0.05) |>
  anti_join(gene_sig, by = c("exposure", "feature_id"))

iso_sum_df <- da_res_exposures |>
  filter(grepl("RNA|Isoforms", exp_name), pvalue < 0.05) |>
  dplyr::count(exposure, exp_name) |>
  mutate(level = if_else(exp_name %in% isoform_assays, "isoform", "gene_protein")) |>
  summarise(
    n_isoform = sum(n[level == "isoform"]),
    n_total    = sum(n),
    .by = exposure
  ) |>
  # Replace raw isoform count with isoform-only count
  left_join(
    isoform_only |> dplyr::count(exposure, name = "n_isoform_only"),
    by = "exposure"
  ) |>
  mutate(
    n_isoform_only = coalesce(n_isoform_only, 0L),
    isoform_frac   = n_isoform_only / n_total
  ) |>
  inner_join(aw_cb, by = c("exposure" = "variable")) |>
  mutate(
    clean_name = tidytext::reorder_within(
      x      = clean_name,
      by     = isoform_frac,
      within = category
    ),
    category = forcats::fct_reorder(
      category, isoform_frac, .fun = median, .desc = TRUE
    )
  )



iso_sum_df |>
  ggplot(aes(x = isoform_frac,
             y = clean_name,
             color = isoform_frac,
             fill  = isoform_frac)) +
  geom_segment(aes(
    x    = 0,
    xend = isoform_frac,
    y    = clean_name
  )) +
  geom_point(shape = 21,
             size   = 2.5,
             stroke = 0.8) +
  tidytext::scale_y_reordered() +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space  = "free",
    strip  = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(iso_sum_df$category)] |>
          scales::alpha(.3)
      )
    )
  ) +
  scale_color_gradientn(
    colors = c("#000004FF", "#B63679FF", "#FCFDBFFF") |> rev()
  ) +
  scale_fill_gradientn(
    colors = c("#000004FF", "#B63679FF", "#FCFDBFFF") |> rev()
  ) +
  scale_x_continuous(labels = scales::percent) +
  geom_vline(xintercept = 0.15, linetype = "dashed", colour = "grey40") +
  labs(
    x     = "Isoform Only %",
    y     = "",
    fill  = "Isoform Only %",
    color = "Isoform Only %"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor        = element_blank(),
    legend.position         = "right",
    legend.title.position   = "top",
    strip.text.y            = element_text(angle = 0, face = "bold.italic")
  )

### Isoform Impact Per Category ---------
iso_sum_df |> 
group_by(category) |> 
  reframe(med_iso_frac=median(isoform_frac))  |> 
  ggplot(aes(
    x=med_iso_frac,
    y=reorder(category,med_iso_frac),
    fill=category
  ))+
  geom_col(width=.9)+
  theme_custom()+
  scale_fill_manual(values = cat_colors)+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(.1,0),
                     labels = scales::percent)+
  labs(
    x="Median Isoform %",
    y=""
  )

## Isoform Impact Per exp_name --------------

# Pair each isoform assay with its RNA counterpart
assay_pairs <- tibble(
  isoform_assay = c("CD16 Monocyte Isoforms", "CD4 T cell Isoforms"),
  rna_assay     = c("CD16 Monocyte RNA",      "CD4 T cell RNA")
)

# Gene-level significant hits, keyed by their paired isoform assay name
gene_sig_by_exp <- assay_pairs |>
  inner_join(
    da_res_exposures |>
      filter(exp_name %in% assay_pairs$rna_assay, pvalue < 0.05) |>
      dplyr::select(exposure, feature_id, rna_assay = exp_name) |>
      distinct(),
    by = "rna_assay"
  ) |>
  dplyr::select(exp_name = isoform_assay, exposure, feature_id)

# Isoform-only per exp_name
isoform_only_by_exp <- da_res_exposures |>
  filter(exp_name %in% assay_pairs$isoform_assay, pvalue < 0.05) |>
  anti_join(gene_sig_by_exp, by = c("exp_name", "exposure", "feature_id"))

iso_sum_df_by_exp <- assay_pairs |>
  pivot_longer(
    everything(),
    names_to  = "level",
    values_to = "exp_name"
  ) |>
  mutate(level = if_else(level == "isoform_assay", "isoform", "gene")) |>
  left_join(
    da_res_exposures |>
      filter(grepl("RNA|Isoforms", exp_name), pvalue < 0.05) |>
      dplyr::count(exposure, exp_name),
    by = "exp_name"
  ) |>
  # Tag each row with its paired isoform assay name before collapsing
  left_join(
    assay_pairs |> rename(exp_name = rna_assay),
    by = "exp_name"
  ) |>
  mutate(
    isoform_assay = if_else(level == "isoform", exp_name, isoform_assay)
  ) |>
  summarise(
    n_isoform = sum(n[level == "isoform"], na.rm = TRUE),
    n_total    = sum(n, na.rm = TRUE),
    .by = c(exposure, isoform_assay)
  ) |>
  rename(exp_name = isoform_assay) |>
  left_join(
    isoform_only_by_exp |>
      dplyr::count(exposure, exp_name, name = "n_isoform_only"),
    by = c("exposure", "exp_name")
  ) |>
  mutate(
    n_isoform_only = coalesce(n_isoform_only, 0L),
    isoform_frac   = n_isoform_only / n_total
  ) |>
  inner_join(aw_cb, by = c("exposure" = "variable")) |>
  mutate(
    clean_name = tidytext::reorder_within(
      x      = clean_name,
      by     = isoform_frac,
      within = interaction(category, exp_name)
    ),
    category = forcats::fct_reorder(
      category, isoform_frac, .fun = median, .desc = TRUE
    )
  )


iso_sum_df_by_exp |>
  ggplot(aes(
    x     = isoform_frac,
    y     = clean_name,
    color = isoform_frac,
    fill  = isoform_frac
  )) +
  geom_segment(aes(
    x    = 0,
    xend = isoform_frac,
    y    = clean_name
  )) +
  geom_point(
    shape  = 21,
    size   = 2.5,
    stroke = 0.8
  ) +
  tidytext::scale_y_reordered() +
  ggh4x::facet_grid2(
    category ~ exp_name,
    scales = "free",
    space  = "free_y",
    strip  = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(iso_sum_df_by_exp$category)] |>
          scales::alpha(.3)
      )
    )
  ) +
  scale_color_gradientn(
    colors = c("#000004FF", "#B63679FF", "#FCFDBFFF") |> rev()
  ) +
  scale_fill_gradientn(
    colors = c("#000004FF", "#B63679FF", "#FCFDBFFF") |> rev()
  ) +
  scale_x_continuous(labels = scales::percent) +
  geom_vline(xintercept = 0.15, linetype = "dashed", colour = "grey40") +
  labs(
    x     = "Isoform Only %",
    y     = "",
    fill  = "Isoform Only %",
    color = "Isoform Only %"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor      = element_blank(),
    legend.position       = "right",
    legend.title.position = "top",
    strip.text.y          = element_text(angle = 0, face = "bold.italic")
  )


### Isoform Impact Per Category × exp_name ---------
iso_sum_df_by_exp |>
  summarise(
    med_iso_frac = median(isoform_frac),
    .by = c(category, exp_name)
  ) |>
  ggplot(aes(
    x    = med_iso_frac,
    y    = reorder(category, med_iso_frac),
    fill = category
  )) +
  geom_col(width = .9) +
  facet_wrap(~ exp_name) +
  theme_custom() +
  scale_fill_manual(values = cat_colors) +
  theme(legend.position = "none") +
  scale_x_continuous(
    expand = c(.1, 0),
    labels = scales::percent
  ) +
  labs(
    x = "Median Isoform %",
    y = ""
  )
## TF Impact --------------
source("./scripts/get_viper.R")

cd4_viper_res <- map(
  unique(da_res_exposures$exposure),~{
    da_res_exposures |> 
      filter(exposure==.x) |> 
      filter(exp_name == "CD4 T cell RNA") |> 
      get_viper() |> 
      mutate(exposure=.x) |> 
      mutate(exp_name="CD4 T cell RNA")
  }) |> 
  bind_rows()

cd16_viper_res <- map(
  unique(da_res_exposures$exposure),~{
    da_res_exposures |> 
      filter(exposure==.x) |> 
      filter(exp_name == "CD16 Monocyte RNA") |> 
      get_viper() |> 
      mutate(exposure=.x) |> 
      mutate(exp_name="CD16 Monocyte RNA")
  }) |> 
  bind_rows()


viper_res <- bind_rows(cd4_viper_res,cd16_viper_res)

cd16_viper_res |> 
  filter(fdr<0.05) |> 
  group_by(tf) |> 
  reframe(
    med_nes = median(nes),
    n_dir = n_distinct(sign(nes)),
    med_logp = median(-log10(fdr)),
    n=n()
  ) |> 
  view()


cd4_viper_res |> 
  filter(fdr<0.05) |> 
  group_by(tf) |> 
  reframe(
    med_nes = median(nes),
    n_dir = n_distinct(sign(nes)),
    med_logp = median(-log10(fdr)),
    n=n()
  ) |> 
  view()


common_tfs <- intersect(
  cd16_viper_res |> 
    filter(fdr<0.05) |> 
    group_by(tf) |> 
    reframe(
      med_nes = median(nes),
      n_dir = n_distinct(sign(nes)),
      med_logp = median(-log10(fdr)),
      n=n()
    ) |> 
    filter(n>20) |> 
    pull(tf),
  cd4_viper_res |> 
    filter(fdr<0.05) |> 
    group_by(tf) |> 
    reframe(
      med_nes = median(nes),
      n_dir = n_distinct(sign(nes)),
      med_logp = median(-log10(fdr)),
      n=n()
    ) |> 
    filter(n>20) |> 
    pull(tf))

cd4_viper_res |> 
  filter(fdr<0.05) |> 
  group_by(tf) |> 
  reframe(
    med_nes = median(nes),
    n_dir = n_distinct(sign(nes)),
    med_logp = median(-log10(fdr)),
    n=n()
  ) |>
  mutate(exp_name = "CD4") |> 
  bind_rows(
    cd16_viper_res |> 
      filter(fdr<0.05) |> 
      group_by(tf) |> 
      reframe(
        med_nes = median(nes),
        n_dir = n_distinct(sign(nes)),
        med_logp = median(-log10(fdr)),
        n=n()
      ) |>
      mutate(exp_name = "CD16") 
  ) |> 
  filter(tf %in% common_tfs) |>
  view()

### Comparison with Lung Function -----------------

cd4_viper_res |> 
  dplyr::rename(exposure_nes = nes) |> 
  inner_join(lung_func_viper_res |> 
               filter(exp_name == "CD4 T cell RNA") |> 
               dplyr::select(tf,lung_func_nes = nes),
             by="tf") |> 
  filter(sign(exposure_nes) != sign(lung_func_nes)) |> 
  inner_join(aw_cb,
             by=c("exposure"="variable")) |> 
  dplyr::count(category) |> 
  view()

cd16_viper_res |> 
  dplyr::rename(exposure_nes = nes) |> 
  inner_join(lung_func_viper_res |> 
               filter(exp_name == "CD16 Monocyte RNA") |> 
               dplyr::select(tf,lung_func_nes = nes),
             by="tf") |> 
  filter(sign(exposure_nes) != sign(lung_func_nes)) |> 
  inner_join(aw_cb,
             by=c("exposure"="variable")) |> 
  dplyr::count(category) |> 
  view()

### TF Plotting ------------
cd4_viper_res |> 
  filter(fdr<0.05) |> 
  group_by(tf) |> 
  reframe(
    med_nes = median(nes),
    n_dir = n_distinct(sign(nes)),
    med_logp = median(-log10(fdr)),
    n=n()
  ) |>
  mutate(exp_name = "CD4 T cell RNA") |> 
  bind_rows(
    cd16_viper_res |> 
      filter(fdr<0.05) |> 
      group_by(tf) |> 
      reframe(
        med_nes = median(nes),
        n_dir = n_distinct(sign(nes)),
        med_logp = median(-log10(fdr)),
        n=n()
      ) |>
      mutate(exp_name = "CD16 Monocyte RNA") 
  ) |> 
  filter(tf %in% common_tfs) |>
  ggplot(aes(
    x=reorder(exp_name,med_nes),
    y=reorder(tf,med_nes),
    fill=med_nes
  ))+
  geom_tile()+
  theme_custom()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="grey85",
    high="#9B6981FF",
    midpoint = 0
  )+
  labs(
    x="",
    y="",
    fill="Median NES"
  )+
  theme(
    axis.text.y = element_text(face="italic")
  )+
  scale_x_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA")
  ))

## CLR + Exposure Association Testing ---------------
decon <- readRDS("./results/da_res/decon.rds")

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
    fis
  )

model_df <- bind_cols(ps, as.data.frame(clr_mat))

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
    fev1_fvc_expom |> 
      extract_results("codebook"),
    by=c("exposure"="variable")
  )

model_df <- model_df |>
  mutate(
    th2_skew = clr_Memory.CD4.T.cell.Th2 -
      ((clr_Memory.CD4.T.cell.Th1 +
          clr_Memory.CD4.T.cell.Th17 
        #clr_Memory.CD4.T.cell.Th1.Th17
      ) / 2)
  ) |> 
  mutate(
    treg_th17 = clr_T.reg - clr_Memory.CD4.T.cell.Th17
  ) |> 
  mutate(
    treg_th2 = clr_T.reg - clr_Memory.CD4.T.cell.Th2
  )

# Then one model per exposure on this single index
th2_skew_res <- map_dfr(num_exp_vars, \(exp_var) {
  lm(
    as.formula(paste("th2_skew ~", exp_var, "+ gli_age + gli_sex + fis")),
    data = model_df
  ) |>
    tidy() |>
    filter(term == exp_var) |>
    mutate(exposure = exp_var)
}) |>
  mutate(fdr = p.adjust(p.value, method = "fdr")) |>
  left_join(aw_cb, by = c("exposure" = "variable"))


treg_th2_res <- map_dfr(num_exp_vars, \(exp_var) {
  lm(
    as.formula(paste("treg_th2 ~", exp_var, "+ gli_age + gli_sex + fis")),
    data = model_df
  ) |>
    tidy() |>
    filter(term == exp_var) |>
    mutate(exposure = exp_var)
}) |>
  mutate(fdr = p.adjust(p.value, method = "fdr")) |>
  left_join(aw_cb, by = c("exposure" = "variable"))


treg_th17_res <- map_dfr(num_exp_vars, \(exp_var) {
  lm(
    as.formula(paste("treg_th17 ~", exp_var, "+ gli_age + gli_sex + fis")),
    data = model_df
  ) |>
    tidy() |>
    filter(term == exp_var) |>
    mutate(exposure = exp_var)
}) |>
  mutate(fdr = p.adjust(p.value, method = "fdr")) |>
  left_join(aw_cb, by = c("exposure" = "variable"))

## Th Bootstrapping --------------
library(rsample)
fit_decon_models <- function(df) {
  map_dfr(
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
            data = df
          ) |>
            tidy() |>
            filter(term == exp_var) |>
            mutate(
              cell     = sub("^clr_", "", clr_col),
              exposure = exp_var
            )
        }
      )
    }
  )
}

set.seed(42)
n_boot <- 100

boot_sign <- rsample::bootstraps(model_df, times = n_boot) |>
  pull(splits) |>
  map_dfr(
    \(split) fit_decon_models(rsample::analysis(split)),
    .id = "boot_id"
  ) |>
  summarise(
    sign_consistency = if (all(is.na(estimate))) {
      NA_real_
    } else {
      max(
        mean(estimate > 0, na.rm = TRUE),
        mean(estimate < 0, na.rm = TRUE)
      )
    },
    dominant_sign = case_when(
      all(is.na(estimate))              ~ NA_character_,
      mean(estimate > 0, na.rm = TRUE) >= 0.5 ~ "positive",
      .default = "negative"
    ),
    .by = c(cell, exposure)
  )

### Add bootstrap res back in ---------------
decon_exp_res <- decon_exp_res |>
  left_join(boot_sign, by = c("cell", "exposure"))


decon_exp_res |> 
  filter(sign_consistency>=0.75) |> 
  dplyr::count(category) |> 
  dplyr::select(category,n_cons=n) |> 
  inner_join(
    aw_cb |> 
      filter(variable %in% num_exp_vars) |> 
      dplyr::count(category) |> 
      dplyr::select(category,n_vars=n),
    by="category"
  ) |> 
  mutate(normed_cons=n_cons/n_vars) |> 
  arrange(desc(normed_cons))


decon_exp_res |> 
  group_by(cell_clean,category) |> 
  reframe(med_est=median(estimate),
          med_logp=median(-log10(p.value)),
          n_sig=length(p.value[p.value<0.1]),
          n_exposures=n_distinct(exposure),
          med_cons=median(sign_consistency)) |> 
  mutate(over_75 = ifelse(med_cons>.75,
                          "Over 75% Consistency",
                          FALSE)) |> 
  ggplot(aes(
    x=reorder(cell_clean,med_est),
    y=reorder(category,med_est),
    fill=med_est
  ))+
  geom_tile()+
  geom_point(aes(size = over_75),
             shape = 22,
             fill = "grey95", 
             color = "black") +
  scale_size_manual(
    values = c("Over 75% Consistency" = 3,
               "FALSE" = NULL))+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0
  )+
  theme_custom()+
  labs(
    x="",
    y="",
    fill="Median Estimate",
    size="Consistency"
  )

library(ggridges)
decon_exp_res |> 
  #filter(p.value<0.1) |> 
  ggplot(aes(
    x=estimate,
    y=fct_reorder(cell_clean,estimate,.fun = median)
  ))+
  geom_density_ridges()+
  xlim(c(-2,2))+
  geom_vline(xintercept = 0,linetype="dashed")

## How many are significant --------------
decon_exp_res |> 
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

## How many are consistent? --------------
decon_exp_res |> 
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

## Th2/Terminal Effector Plotting ------------
decon_exp_res |> 
  filter(sign_consistency>=0.75) |>
  filter(p.value<0.1) |> 
  filter(grepl("Th2|effector",cell)) |> 
  #mutate(clean_name=tidytext::reorder_within(clean_name,by = estimate,within = cell_clean)) |> 
  ggplot(aes(
    x=reorder(cell_clean,estimate),
    y=reorder(clean_name,estimate),
    fill=estimate
  ))+
  geom_tile()+
  theme_custom()+
  #tidytext::scale_y_reordered()+
  #facet_grid(~cell_clean~.,scales="free",space="free")+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0
  )+
  #scale_x_continuous(expand = c(.1,0))+
  labs(
    x="",
    y="",
    fill="Estimate"
  )
  # theme(
  #   axis.text.x = element_text(angle=0,hjust=.6)
  # )
  

### Th - Exposure Plotting --------
library(patchwork)

decon_exp_res |>
  #filter(!grepl("Serum Essential|Urine Essential",category)) |> 
  mutate(z = sign(estimate) * qnorm(1 - p.value / 2)) |>
  summarise(
    stouffer_z = sum(z) / sqrt(n()),
    stouffer_p = 2 * pnorm(-abs(stouffer_z)),
    n_tests    = n(),
    .by = c(cell_clean, category)
  ) |>
  mutate(
    category = reorder(category, stouffer_z, FUN = median),
    cell_clean = reorder(cell_clean, stouffer_z, FUN = median)
  ) |>
  (\(d) {
    
    main <- d |>
      ggplot(aes(
        x    = category,
        y    = cell_clean,
        fill = stouffer_z,
        size = -log10(stouffer_p)
      )) +
      geom_point(shape = 21) +
      theme_custom() +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits=c(-4,4),
        oob=scales::squish
      ) +
      labs(
        x    = "",
        y    = "",
        fill = "Stouffer Z",
        size=expression("-Log"[10]*"P")
      ) +
      theme(axis.text.x = element_text(angle = 90, vjust = .6))
    
    top <- d |>
      summarise(med_z = median(stouffer_z), .by = category) |>
      mutate(category = factor(category, levels = levels(d$category))) |>
      ggplot(aes(
        x    = category,
        y    = abs(med_z),
        fill = med_z
      )) +
      geom_col() +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits=c(-4,4),
        oob=scales::squish
      ) +
      theme_custom() +
      labs(x = "", 
           y = "|Median Z|",
           fill="Median Z") +
      theme(
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "none"
      )
    
    right <- d |>
      summarise(med_z = median(stouffer_z), .by = cell_clean) |>
      mutate(cell_clean = factor(cell_clean, levels = levels(d$cell_clean))) |>
      ggplot(aes(
        x    = abs(med_z),
        y    = cell_clean,
        fill = med_z
      )) +
      geom_col() +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits=c(-4,4),
        oob=scales::squish
      ) +
      theme_custom() +
      labs(x = "|Median Z|",
           y = "",
           fill = "Median Z") +
      theme(
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        legend.position = "none"
      )
    
    patchwork::wrap_plots(
      top, patchwork::plot_spacer(),
      main, right,
      ncol    = 2,
      widths  = c(4, 1),
      heights = c(1.5, 4)
    )+plot_layout(guides = "collect")
    
  })()

### Category Th2:Th1/17 and Treg Ratio ----------
th2_skew_res |> 
  group_by(category) |> 
  reframe(
    med_est= median(estimate),
    med_logp=median(-log10(p.value)),
    n=n()
  ) |> 
  mutate(comparison="Th2:Th1/17") |> 
  bind_rows(
    treg_th2_res |> 
      group_by(category) |> 
      reframe(
        med_est= median(estimate),
        med_logp=median(-log10(p.value)),
        n=n()
      ) |> 
      mutate(comparison="Treg:Th2")
  ) |> 
  bind_rows(
    treg_th17_res |> 
      group_by(category) |> 
      reframe(
        med_est= median(estimate),
        med_logp=median(-log10(p.value)),
        n=n()
      ) |> 
      mutate(comparison="Treg:Th17")
  ) |> 
  ggplot(aes(
    x=reorder(comparison,med_est),
    y=reorder(category,med_est),
    fill=med_est,
    size=med_logp
  ))+
  geom_point(shape=21)+
  theme_custom()+
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "white",
    high     = "#9B6981FF",
    midpoint = 0)+
  labs(
    x="",
    y="",
    size=expression("Median -Log"[10]*"P"),
    fill="Median Estimate"
  )


### Exposures Th2:Th1/17 and Treg Ratio ----------
# Build the combined summary with clean_name and category joined in
combined_skew_res <- list(
  "Th2:Th1/17" = th2_skew_res,
  "Treg:Th2"   = treg_th2_res,
  "Treg:Th17"  = treg_th17_res
) |>
  purrr::imap(\(df, comp) df |> mutate(comparison = comp)) |>
  purrr::list_rbind() |>
  summarise(
    med_est  = median(estimate),
    med_logp = median(-log10(p.value)),
    .by = c(clean_name, category, comparison)
  ) |>
  mutate(
    category  = forcats::fct_reorder(category, med_est, .fun = median, .desc = TRUE),
    clean_name = tidytext::reorder_within(clean_name, med_est, within = category)
  )

combined_skew_res |>
  ggplot(aes(
    x    = reorder(comparison, med_est),
    y    = clean_name,
    fill = med_est,
    size = med_logp
  )) +
  geom_point(shape = 21) +
  tidytext::scale_y_reordered() +
  ggh4x::facet_grid2(
    category ~ .,
    scales = "free",
    space  = "free",
    strip  = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = cat_colors[levels(combined_skew_res$category)] |>
          scales::alpha(.3)
      )
    )
  ) +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "white",
    high     = "#9B6981FF",
    midpoint = 0,
    limits=c(-2,2),
    oob =scales::squish
  ) +
  theme_custom() +
  theme(
    strip.text.y          = element_text(angle = 0, face = "bold.italic"),
    panel.grid.minor      = element_blank(),
    legend.position       = "right",
    legend.title.position = "top"
  ) +
  labs(
    x    = "",
    y    = "",
    size = expression("-Log"[10]*"P"),
    fill = "Estimate"
  )


## Top Exposure Driven Pathways ------------

### Old - CD16 ----------
exp_gsea_res_sig_clean |>
  filter(size < 200) |>
  filter(exp_name == "CD16 Monocyte RNA") |>
  group_by(term_name) |>
  reframe(
    med_nes    = median(NES),
    med_logp   = median(-log10(padj)),
    n_exposures = n_distinct(exposure),
    n          = n(),
    n_dir      = n_distinct(sign(NES)),
    n_up       = length(exposure[NES > 0]),
    n_down     = length(exposure[NES < 0]),
    categories = paste(unique(category), collapse = ","),
    n_cat      = n_distinct(category)
  ) |>
  arrange(desc(n_exposures)) |>
  slice_head(n = 10) |>
  arrange(desc(med_nes)) |> 
  mutate(term_name = factor(term_name)) |> 
  (\(df) {
    
    shared_y <- scale_y_discrete(
      limits = df$term_name[order(df$med_nes)])
    
    p1 <- df |>
      ggplot(aes(x = med_nes, 
                 y = term_name,
                 fill = med_nes)) +
      geom_col() +
      shared_y +
      scale_fill_gradient2(
        low="#5399b0",
        mid="grey85",
        high="#9B6981FF",
        midpoint = 0,
        limits=c(-2,2),
        oob=scales::squish
      )+
      labs(x = "Median NES",
           y = "",
           fill = "Median NES",
           title = expression(bold("CD16")^bold("+") * bold(" Monocyte RNA"))) +
      theme_custom()
    
    p2 <- df |>
      ggplot(aes(x = 1, 
                 y = term_name,
                 size =n_exposures)) +
      geom_point(shape=21,fill="grey38") +
      shared_y +
      labs(x = "",
           y = "",
           size="No. of Exposure") +
      theme_custom() +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    
    p1 + p2 + plot_layout(nrow = 1,
                          widths = c(2,1),
                          guides = "collect")
    
  })()

### Old - All Omics ----------
exp_gsea_res_sig_clean |>
  filter(size < 200) |>
  group_by(exp_name, term_name) |>
  reframe(
    med_nes     = median(NES),
    med_logp    = median(-log10(padj)),
    n_exposures = n_distinct(exposure),
    n           = n(),
    n_dir       = n_distinct(sign(NES)),
    n_up        = length(exposure[NES > 0]),
    n_down      = length(exposure[NES < 0]),
    pct_up      = n_up / (n_up + n_down),
    categories  = paste(unique(category), collapse = ","),
    n_cat       = n_distinct(category)
  ) |>
  group_by(exp_name) |>
  slice_max(n_exposures, n = 8, with_ties = F) |>
  ungroup() |>
  mutate(
    term_name = tidytext::reorder_within(term_name, med_nes, exp_name),
    is_mixed  = n_dir > 1
  ) |>
  ggplot(aes(
    x    = med_nes,
    y    = term_name,
    fill = med_nes
  )) +
  geom_col() +
  geom_point(
    aes(
      x      = max(med_nes) * 1.25,
      size   = n_exposures,
      colour = pct_up,
      shape  = is_mixed
    )
  ) +
  scale_shape_manual(
    values = c("FALSE" = 21, "TRUE" = 23),  # circle vs diamond
    labels = c("FALSE" = "Consistent", "TRUE" = "Mixed direction")
  ) +
  scale_colour_gradient2(
    low      = "#5399b0",
    mid      = "grey85",
    high     = "#9B6981FF",
    midpoint = 0.5,
    limits   = c(0, 1),
    labels   = scales::percent
  ) +
  tidytext::scale_y_reordered() +
  scale_fill_gradient2(
    low      = "#5399b0",
    mid      = "grey85",
    high     = "#9B6981FF",
    midpoint = 0,
    limits   = c(-2, 2),
    oob      = scales::squish
  ) +
  scale_x_continuous(expand = c(.1,0))+
  facet_wrap(
    ~ exp_name,
    scales   = "free_y",
    #space    = "free_y",
    nrow     = 1,
    labeller = as_labeller(c(
      "Adduct Burden"     = "bold('Adduct Burden')",
      "CD16 Monocyte RNA" = "bold('CD16')^bold('+')~bold('Monocyte RNA')",
      "CD4 T cell RNA"    = "bold('CD4')^bold('+')~bold('T cell RNA')",
      "Proteomics"        = "bold('Proteomics')"
    ), default = label_parsed)
  )+
  labs(
    x      = "Median NES",
    y      = "",
    fill   = "Median NES",
    size   = "No. of Exposures",
    colour = "% Exposures Up",
    shape  = "Direction"
  ) +
  theme_custom()+
  theme(
    legend.position = "bottom",
    legend.title.position = "top"
  )

### Current GSEA Plot --------
source("./scripts/04_1_exposure_omics_assoc_gsea_plot.R")
(p1 / p2 / p3 / p4) +
  plot_layout(guides = "collect", heights = c(3, 1, 1, 1))&
  theme(legend.position = "bottom")

### Exposures Driving GSEA ----------
# bp1&theme(legend.position = "none")
# bp2+
#   scale_y_break(c(15, 110),ticklabels = c(0,15,110,115)) +
#   theme_bw() +
#   theme_custom()+
#   theme(legend.position = "none",
#         axis.text.x           = element_blank(),
#         axis.ticks.x          = element_blank(),
#         axis.line.x           = element_blank())
# bp3+
#   scale_y_break(c(15, 110),ticklabels = c(0,15,110,115)) +
#   theme_bw() +
#   theme_custom()+
#   theme(legend.position = "none",
#         axis.text.x           = element_blank(),
#         axis.ticks.x          = element_blank(),
#         axis.line.x           = element_blank())
# bp4+
#   scale_y_break(c(15, 110),ticklabels = c(0,15,110,115)) +
#   theme_bw() +
#   theme_custom()+
#   theme(legend.position = "none",
#         axis.text.x           = element_blank(),
#         axis.ticks.x          = element_blank(),
#         axis.line.x           = element_blank())

bp1&guides(fill="none")
bp2&guides(fill="none")
bp3&guides(fill="none")
bp4&guides(fill="none")

# now just one to grab the legend
bp1&theme(legend.position = "bottom")

## Consistent Exposure-omics pathways -----------
exp_gsea_res_sig_clean |>
  group_by(term_name) |>
  reframe(
    n_dir      = n_distinct(sign(NES)),
    n_exposure = n_distinct(exposure),
    med_nes    = median(NES),
    n_omics    = n_distinct(exp_name)
  ) |>
  filter(n_dir == 1) |>
  arrange(med_nes) |>
  mutate(
    rank = row_number(),
    label = case_when(
      term_name %in% c(
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
      ) ~ term_name,
      .default = NA
    )
  ) |>
  (\(df) {
    label_df <- df |>
      filter(!is.na(label)) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(3, 12, length.out = n()))
    
    df |>
      mutate(exp_name = "Exposure GSEA Terms") |>
      ggplot(aes(x = rank, y = exp_name, fill = med_nes)) +
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
        aes(x = x_pos + 0.6, y = y_label, label = label),
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
      theme(
        axis.text.x  = element_blank(),
        axis.line.y  = element_blank(),
        axis.line.x  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title.position = "top",
        plot.margin  = margin(40, 200, 5, 5)
      ) +
      labs(x = NULL, 
           y = NULL, 
           fill = "Median NES")
  })()

## Regulons per exposure category -----------

list(
  cd16 = cd16_viper_res,
  cd4  = cd4_viper_res
) |>
  imap(\(df, nm) df |> mutate(source = nm)) |>
  list_rbind() |>
  left_join(aw_cb, by = c("exposure" = "variable")) |>
  filter(fdr < 0.05) |>
  group_by(source, category) |>
  reframe(
    n_exposures = n_distinct(exposure),
    n_tf        = n_distinct(tf)
  ) |>
  mutate(
    source = factor(source, levels = c("cd16", "cd4")),
    category_ordered = tidytext::reorder_within(category, n_tf, source)
  ) |>
  ggplot(aes(
    x    = n_tf,
    y    = category_ordered,
    fill = category
  )) +
  geom_col(width = 0.9) +
  tidytext::scale_y_reordered() +
  scale_fill_manual(values = cat_colors) +
  facet_wrap(
    ~ source,
    scales = "free_y",
    labeller = as_labeller(c(
      cd16 = "CD16\u207a Monocyte RNA",
      cd4  = "CD4\u207a T cell RNA"
    ))
  ) +
  theme_custom() +
  guides(fill = "none") +
  labs(x = "No. of Regulons", y = NULL)



### Regulon Bubble plot --------
list(
  cd16 = cd16_viper_res,
  cd4  = cd4_viper_res
) |>
  imap(\(df, nm) df |> mutate(source = nm)) |>
  list_rbind() |>
  left_join(aw_cb, by = c("exposure" = "variable")) |>
  filter(fdr < 0.05) |>
  group_by(source, category) |>
  reframe(
    n_exposures = n_distinct(exposure),
    n_tf        = n_distinct(tf)
  ) |>
  pivot_wider(names_from = source,values_from = n_tf) |> 
  ggplot(aes(
    x=cd16,
    y=cd4,
    fill=category
  ))+
  geom_point(shape=21,size=6)+
  scale_fill_manual(values=cat_colors)+
  theme_custom()+
  labs(
    x=expression("No. of CD16"^"+" ~ "Monocyte Regulons"),
    y=expression("No. of CD4"^"+" ~ "T cell Regulons"),
    fill="Category"
  )
  

## No. of Pubmed hits ----------------

go_pct_res <- map(unique(da_res_exposures$exposure), ~{
  
  cd4_genes <- da_res_exposures |>
    filter(pvalue < 0.05) |>
    filter(exposure == .x) |>
    filter(grepl("CD4", exp_name)) |>
    pull(feature_id) |>
    unique()
  
  cd4_go_tbl <- go_term_data$mapping |>
    filter(gene_symbol %in% cd4_genes) |>
    dplyr::count(gene_symbol, name = "n_go")
  
  cd4_go_tbl <- cd4_go_tbl |>
    bind_rows(
      data.frame(
        gene_symbol = cd4_genes[!cd4_genes %in% cd4_go_tbl$gene_symbol],
        n_go        = 0
      )
    ) |>
    mutate(over_5 = n_go > 5)
  
  cd4_pct <- cd4_go_tbl |>
    dplyr::count(over_5) |>
    mutate(
      total        = nrow(cd4_go_tbl),
      pct_below_5  = n / total
    ) |>
    filter(over_5 == FALSE) |>
    dplyr::select(cd4_pct_below_5 = pct_below_5)
  
  
  cd16_genes <- da_res_exposures |>
    filter(pvalue < 0.05) |>
    filter(exposure == .x) |>
    filter(grepl("CD16", exp_name)) |>
    pull(feature_id) |>
    unique()
  
  cd16_go_tbl <- go_term_data$mapping |>
    filter(gene_symbol %in% cd16_genes) |>
    dplyr::count(gene_symbol, name = "n_go")
  
  cd16_go_tbl <- cd16_go_tbl |>
    bind_rows(
      data.frame(
        gene_symbol = cd16_genes[!cd16_genes %in% cd16_go_tbl$gene_symbol],
        n_go        = 0
      )
    ) |>
    mutate(over_5 = n_go > 5)
  
  cd16_pct <- cd16_go_tbl |>
    dplyr::count(over_5) |>
    mutate(
      total       = nrow(cd16_go_tbl),
      pct_below_5 = n / total
    ) |>
    filter(over_5 == FALSE) |>
    dplyr::select(cd16_pct_below_5 = pct_below_5)
  
  
  data.frame(
    exposure         = .x,
    cd4_pct_below_5  = if (nrow(cd4_pct) > 0)  cd4_pct$cd4_pct_below_5   else NA_real_,
    cd16_pct_below_5 = if (nrow(cd16_pct) > 0) cd16_pct$cd16_pct_below_5 else NA_real_
  )
  
}) |>
  bind_rows() |>
  left_join(
    aw_cb |> filter(variable %in% num_exp_vars),
    by = c("exposure" = "variable")
  )

### Bubble Plot -----------
go_pct_res |> 
  group_by(category) |> 
  reframe(
    med_cd4_pct = median(cd4_pct_below_5),
    med_cd16_pct = median(cd16_pct_below_5)
  ) |> 
  ggplot(aes(
    x=med_cd16_pct,
    y=med_cd4_pct,
    fill=category
  ))+
  geom_point(shape=21,size=6)+
  scale_fill_manual(values=cat_colors)+
  theme_custom()+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  labs(
    x=expression("Median CD16"^"+" ~ "Monocyte"[Understudied]),
    y=expression("Median CD4"^"+" ~ "T cell"[Understudied]),
    fill="Category"
  )
  

## Save Data ----------------------
# saveRDS(da_res_exposures,file = "./results/exposure_omics/da_res_exposures.rds")

saveRDS(decon_exp_res,file="./results/exposure_omics/decon_exp_res.rds")