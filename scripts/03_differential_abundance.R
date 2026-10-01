## Load Libraries ------------------------
library(tidyverse)
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
cat_colors <- scales::alpha(get_palette("aaas",13), 0.5) |> 
  (\(colors){
    names(colors) <- aw_cb |> 
      filter(variable %in% exp_vars) |> 
      pull(category) |> 
      unique()
    colors
  })()

msig <- msigdbr::msigdbr(species = "human",collection = "C5")

## Normalization Function --------------
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

## Run Differential Exp./Abd. -----------
base_formula <- ~ pftfev1fvc_actual + gli_age + gli_sex + fis

da_res_total <- bind_rows(
  run_deseq(mae = fev1_fvc_expom,
            experiment_name = "CD4 T cell RNA",
            fm = base_formula),
  run_deseq(mae = fev1_fvc_expom,
            experiment_name = "CD16 Monocyte RNA",
            fm = base_formula),
  run_deseq(mae = fev1_fvc_expom,
            experiment_name = "CD4 T cell Isoforms",
            fm = base_formula),
  run_deseq(mae = fev1_fvc_expom,
            experiment_name = "CD16 Monocyte Isoforms",
            fm = base_formula),
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



### No. of DEGs ----------------------
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
    title = bquote(bold("FEV")[1] ~ bold("/FVC")),
    subtitle = "",
    caption = "Nominal P < 0.05"
  )

### Volcano Plot ------------

da_res_total |>
  filter(!is.na(feature_id)) |> 
  filter(!grepl("^ENSG",feature_id)) |> 
  mutate(
    sig = pvalue < 0.05,
    direction = case_when(
      pvalue < 0.05 & logfc > 0 ~ "Up",
      pvalue < 0.05 & logfc < 0 ~ "Down",
      .default = "ns"
    ) |> factor(levels = c("Up", "ns", "Down"))
  ) |>
  ggplot(aes(
    x     = logfc,
    y     = -log10(pvalue),
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
      group_by(exp_name) |> 
      slice_min(logfc, n = 5),
    aes(label = feature_id),
    color              = "#5399b0",
    fill               = "white",
    size               = 3,
    fontface           = "italic",
    box.padding        = 0.5,
    max.overlaps = 20,
    min.segment.length = 0,
    seed               = 42
  ) +
  ggrepel::geom_label_repel(
    data = \(d) d |>
      filter(direction == "Up") |>
      group_by(exp_name) |> 
      slice_min(logfc, n = 5),
    aes(label = feature_id),
    color              = "#9B6981FF",
    fill               = "white",
    size               = 3,
    fontface           = "italic",
    box.padding        = 0.5,
    max.overlaps = 20,
    min.segment.length = 0,
    seed               = 123
  ) +
  scale_color_manual(values = c(
    "Up"   = "#9B6981FF",
    "ns"   = "grey75",
    "Down" = "#5399b0"
  )) +
  labs(
    x     = expression("Log"[2]*"FC"),
    y     = expression("-Log"[10]*"P"),
    color = "Direction"
  ) +
  facet_wrap(~exp_name,space="free_y",scales="free_x",ncol = 1)+
  theme_custom()+
  theme(
    axis.text.x=element_text(angle=0,hjust=.6)
  )

## Protein v. Adduct -----------
da_res_total |> 
  filter(exp_name %in% c(
    "Proteomics"
  )) |> 
  dplyr::select(
    "protein_logfc"=logfc,
    feature_id
  ) |> 
  # remove extreme values
  #filter(protein_logfc>-4,protein_logfc<4) |> 
  inner_join(
    da_res_total |> 
      filter(exp_name %in% c(
        "Adduct Burden"
      )) |> 
      dplyr::select(
        "adduct_logfc"=logfc,
        feature_id
      ),
    by="feature_id"
  ) |> 
  ggplot(aes(
    x=protein_logfc,
    y=adduct_logfc
  ))+
  geom_point()+
  theme_custom()+
  geom_point(color="grey35")+
  stat_smooth(method = "lm",
              color="#B63679FF")+
  stat_cor(p.accuracy = 0.001)+
  labs(
    x=expression("Protein Log"[2]*"FC"),
    y=expression("Adduct Log"[2]*"FC")
  )+
  theme(
    axis.text.x=element_text(angle=0,hjust=0.5)
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
    #padj < 0.1
  )


gsea_res_sig <- gsea_res |> 
  filter(padj < 0.1)


## Go Enrichment --------------
da_lst <- da_res_total |> 
  filter(grepl("Adductomics|Isoform",exp_name)) |> 
  filter(pvalue<0.05) |> 
  mutate(direction=ifelse(
    logfc>0,
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

### Enrichment Plotting -----------
gsea_res_sig |> 
  mutate(direction=ifelse(NES>0,"Up","Down")) |> 
  group_by(group,direction) |> 
  arrange(desc(-log10(pval))) |> 
  slice_head(n=10) |> 
  ggplot(aes(
    x=NES,
    y=reorder(str_trunc(term_name,width = 40),NES),
    fill=NES
  ))+
  geom_col()+
  theme_custom()+
  facet_grid(group~.,scales="free_y",space="free_y")+
  scale_fill_gradient2(
    low="#5399b0",
    mid="grey85",
    high="#9B6981FF",
    midpoint = 0
  )+
  theme(
    axis.text.x = element_text(angle=0,hjust=.5)
  )+
  labs(
    x="NES",
    y="",
    fill="NES"
  )


enr_res_sig |> 
  separate(group,into = c("exp_name","direction"),sep = "_",remove = F) |> 
  mutate(group=gsub("_","-",group)) |> 
  group_by(exp_name,direction) |> 
  arrange(desc(-log10(p_adjust))) |> 
  slice_head(n=10) |> 
  ggplot(aes(
    x=-log10(p_adjust),
    y=reorder(str_trunc(term_name,width=40),-log10(p_adjust)),
    fill=-log10(p_adjust)
  ))+
  geom_col()+
  theme_custom()+
  facet_grid(group~.,scales="free_y",space="free_y")+
  scale_fill_gradient(
    high = "#360033", 
    low = "#FFFDE4"
  )+
  theme(
    axis.text.x = element_text(angle=0,hjust=.5)
  )+
  labs(
    x=expression("-Log"[10]*"P"),
    y="",
    fill=expression("-Log"[10]*"P")
  )

### Focused GSEA ----------
gsea_res_sig |>
  filter(group == "CD16 Monocyte RNA") |> 
  mutate(direction=ifelse(NES>0,"Up","Down")) |> 
  group_by(group,direction) |> 
  arrange(desc(-log10(pval))) |> 
  slice_head(n=10) |>
  
  bind_rows(
    gsea_res_sig |>
      filter(group == "Proteomics") |> 
      mutate(direction=ifelse(NES>0,"Up","Down")) |> 
      group_by(group,direction) |> 
      arrange(desc(-log10(pval))) |> 
      slice_head(n=10)
  ) |> 
  
  bind_rows(
    gsea_res_sig |>
      filter(group == "CD4 T cell RNA") |> 
      mutate(direction=ifelse(NES>0,"Up","Down")) |> 
      filter(grepl("mito|tricarboxylic|alpha-beta|chromatin|nucleosome|integrin|adhesion|ubiquitin|SUMO|catbolic",term_name)) 
    
  ) |> 
  ggplot(aes(
    x=NES,
    y=reorder(str_trunc(term_name,width = 40),NES),
    fill=NES
  ))+
  geom_col()+
  theme_custom()+
  facet_grid(group~.,scales="free_y",space="free_y")+
  scale_fill_gradient2(
    low="#5399b0",
    mid="grey85",
    high="#9B6981FF",
    midpoint = 0
  )+
  theme(
    axis.text.x = element_text(angle=0,hjust=.5)
  )+
  labs(
    x="NES",
    y="",
    fill="NES"
  )

### Enr Tile Plots ---------------

#### GSEA ----------------
gsea_res_sig |>
  filter(group == "CD16 Monocyte RNA") |> 
  mutate(direction=ifelse(NES>0,"Up","Down")) |> 
  group_by(group,direction) |> 
  arrange(desc(-log10(pval))) |> 
  slice_head(n=10) |>
  
  bind_rows(
    gsea_res_sig |>
      filter(group == "Proteomics") |> 
      mutate(direction=ifelse(NES>0,"Up","Down")) |> 
      group_by(group,direction) |> 
      arrange(desc(-log10(pval))) |> 
      slice_head(n=10)
  ) |> 
  
  bind_rows(
    gsea_res_sig |>
      filter(group == "CD4 T cell RNA") |> 
      mutate(direction=ifelse(NES>0,"Up","Down")) |> 
      filter(grepl("mito|tricarboxylic|alpha-beta|chromatin|nucleosome|integrin|adhesion|ubiquitin|SUMO|catbolic",
                   term_name)) 
    
  ) |> 
  mutate(term_label = str_trunc(term_name, width = 40)) |>
  mutate(term_label = reorder(term_label, NES)) |> 
  ggplot(aes(
    x=group,
    y=term_label,
    fill=NES
  ))+
  geom_tile()+
  theme_custom()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="grey85",
    high="#9B6981FF",
    midpoint = 0
  )+
  theme(
    axis.text.x = element_text(angle=90,vjust=.6)
  )+
  labs(
    x="NES",
    y="",
    fill="NES"
  )


#### Clean Tile Plots - GSEA Style -------------
# ggplot solution ----
highlight_terms <- c(
  "positive regulation of interleukin-2 production",
  "positive regulation of inflammatory response",
  "defense response to virus",
  
  "immunoglobulin mediated immune response",
  "antigen binding",
  "immunoglobulin complex",
  
  "alpha-beta T cell receptor complex",
  "mitochondrial small ribosomal subunit",
  "integrin-mediated signaling pathway",
  "chromatin remodeling",
  "ubiquitin-dependent protein catabolic process",
  "cellular response to hypoxia",
  "nucleosome assembly"
  
)


p1 <- gsea_res_sig |>
  arrange(desc(NES)) |>
  mutate(rank = rank(NES)) |>
  (\(df) {
    label_df <- df |>
      filter(group == "CD4 T cell RNA", term_name %in% highlight_terms) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(3,10, length.out = n()))
    
    df |>
      mutate(
        fill_NES = if_else(group == "CD4 T cell RNA", NES, NA_real_),
        exp_name = "CD4 T cell RNA"
      ) |>
      ggplot(aes(x = rank, 
                 y = exp_name, 
                 fill = fill_NES)) +
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
        data        = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes = FALSE,
        hjust       = 0,
        size        = 4.5,
        color       = "grey20",
        fill        = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits = c(-2,2),
        oob=scales::squish,
        na.value = "grey92"
      ) +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "Proteomics" ~ "Proteomics"
      ))+
      coord_cartesian(clip = "off") +
      theme_custom() +
      #theme_void()+
      theme(
        axis.text.x  = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(40, 150, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "NES")
  })()

p2 <- gsea_res_sig |>
  arrange(desc(NES)) |>
  mutate(rank = rank(NES)) |>
  (\(df) {
    label_df <- df |>
      filter(group == "CD16 Monocyte RNA", term_name %in% highlight_terms) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(2,5, length.out = n()))
    
    df |>
      mutate(
        fill_NES = if_else(group == "CD16 Monocyte RNA", NES, NA_real_),
        exp_name = "CD16 Monocyte RNA"
      ) |>
      ggplot(aes(x = rank, 
                 y = exp_name, 
                 fill = fill_NES)) +
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
        data        = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes = FALSE,
        hjust       = 0,
        size        = 4.5,
        color       = "grey20",
        fill        = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits = c(-2,2),
        oob=scales::squish,
        na.value = "grey92"
      ) +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "Proteomics" ~ "Proteomics"
      ))+
      coord_cartesian(clip = "off") +
      theme_custom() +
      #theme_void()+
      theme(
        axis.text.x  = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(40, 150, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "NES")
  })()



p3 <- gsea_res_sig |>
  arrange(desc(NES)) |>
  mutate(rank = rank(NES)) |>
  (\(df) {
    label_df <- df |>
      filter(group == "Proteomics", term_name %in% highlight_terms) |>
      mutate(x_pos = rank) |>
      arrange(x_pos) |>
      mutate(y_label = seq(2,5, length.out = n()))
    
    df |>
      mutate(
        fill_NES = if_else(group == "Proteomics", NES, NA_real_),
        exp_name = "Proteomics"
      ) |>
      ggplot(aes(x = rank, 
                 y = exp_name, 
                 fill = fill_NES)) +
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
        data        = label_df,
        aes(x = x_pos + 0.6, y = y_label, label = term_name),
        inherit.aes = FALSE,
        hjust       = 0,
        size        = 4.5,
        color       = "grey20",
        fill        = "white",
        label.padding = unit(0.1, "lines")
      ) +
      scale_fill_gradient2(
        low      = "#5399b0",
        mid      = "white",
        high     = "#9B6981FF",
        midpoint = 0,
        limits = c(-2,2),
        oob=scales::squish,
        na.value = "grey92"
      ) +
      scale_y_discrete(labels=c(
        "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
        "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA"),
        "Proteomics" ~ "Proteomics"
      ))+
      coord_cartesian(clip = "off") +
      theme_custom() +
      #theme_void()+
      theme(
        axis.text.x  = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(40, 150, 5, 5)
      ) +
      labs(x = NULL, y = NULL, fill = "NES")
  })()

library(patchwork)
(p1/p2/p3)+plot_layout(guides = "collect",heights = c(2,1,1))


#### GO ------

enr_res_sig |>
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
    axis.text.x = element_text(angle = 90, vjust = .6)
  ) +
  labs(
    x = "",
    y = "",
    fill = expression("-Log"[10]*"P")
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

rna_assays <- c("CD4 T cell RNA", "CD16 Monocyte RNA",
                "CD4 T cell Isoforms", "CD16 Monocyte Isoforms")

expom_normed <- rna_assays |>
  purrr::reduce(normalize_rna_assay, .init = fev1_fvc_expom)

cs_mat <- tidyexposomics:::.update_assay_colData(
  expom_normed,
  "CD4 T cell RNA") |> 
  assay("counts")

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
    fis,
    .sample
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

## Viper TF ----------------
source("./scripts/get_viper.R")
lung_func_viper_res <- da_res_total |> 
  filter(exp_name =="CD4 T cell RNA") |> 
  get_viper() |> 
  mutate(exp_name = "CD4 T cell RNA") |> 
  bind_rows(
    da_res_total |> 
      filter(exp_name =="CD16 Monocyte RNA") |> 
      get_viper() |> 
      mutate(exp_name = "CD16 Monocyte RNA")
  ) |> 
  filter(fdr<0.05)


lung_func_viper_res |> 
  ggplot(aes(
    y=reorder(exp_name,nes),
    x=reorder(tf,nes),
    fill=nes
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
    fill="NES"
  )+
  theme(
    axis.text.x = element_text(face="italic"),
    # legend.position = "bottom",
    # legend.title.position = "top"
  )+
  scale_y_discrete(labels=c(
    "CD4 T cell RNA"       = expression(CD4^"+" ~ "T cell RNA"),
    "CD16 Monocyte RNA"  = expression(CD16^"+" ~ "Monocyte RNA")
  ))

dorothea <- dorothea::dorothea_hs |> 
  filter(confidence %in% c("A","B","C","D"))

### GO ----------

lung_func_targets <- dorothea |> 
  filter(tf %in% lung_func_viper_res$tf) |> 
  inner_join(lung_func_viper_res,
             by=c("tf")) |> 
  mutate(direction = ifelse(sign(nes) > 0,
                            "Up",
                            "Down")) |> 
  mutate(combo = paste(exp_name,tf,direction,sep="-")) |> 
  (\(df) split(df$target, df$combo))()

# Run enrichment with matched universes
lung_func_enr_res <- future_imap(
  lung_func_targets,
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

lung_func_enr_res_sig <- lung_func_enr_res |>
  filter(
    p_adjust < .1,
    n_with_sel > 5,
    N_with < 1000
  )|> 
  separate(group,sep="-",
           into = c("exp_name","tf","direction"),
           remove = F) 


group_labels <- c(
  "CD16 Monocyte RNA-Down" = "bold(CD16^+~Monocyte~Down)",
  "CD16 Monocyte RNA-Up"   = "bold(CD16^+~Monocyte~Up)",
  "CD4 T cell RNA-Down"    = "bold(CD4^+~T~Cell~Down)",
  "CD4 T cell RNA-Up"      = "bold(CD4^+~T~Cell~Up)"
)

lung_func_enr_res_sig |> 
  group_by(exp_name,direction,term_name) |> 
  reframe(
    n=n(),
    med_logp=median(-log10(p_adjust+1e-20))
  ) |> 
  group_by(exp_name,direction) |> 
  arrange(desc(n),desc(med_logp)) |> 
  slice_head(n=10) |> 
  ungroup() |> 
  mutate(combo=paste(exp_name,direction,sep="-")) |> 
  mutate(term_name = tidytext::reorder_within(
    str_trunc(term_name,width=40),
    by = med_logp,
    within = combo)) |> 
  ggplot(aes(
    x=med_logp,
    y=term_name,
    fill=med_logp
  ))+
  geom_col()+
  theme_custom()+
  scale_fill_gradient(
    high = "#360033",
    low  = "#FFFDE4"
  ) +
  tidytext::scale_y_reordered()+
  facet_wrap(
    ~ tf,
    nrow = 2,
    labeller = as_labeller(
      c(
        "CD16 Monocyte RNA-Down" = "bold(CD16^'+'~Monocyte~RNA~Down)",
        "CD16 Monocyte RNA-Up"   = "bold(CD16^'+'~Monocyte~RNA~Up)",
        "CD4 T cell RNA-Down"    = "bold(CD4^'+'~T~cell~RNA~Down)",
        "CD4 T cell RNA-Up"      = "bold(CD4^'+'~T~cell~RNA~Up)"
      ),
      label_parsed
    ),
    scales = "free_y"
  )+
  labs(
    x=expression("Median -Log"[10]*"P"),
    y="",
    fill=expression("Median -Log"[10]*"P")
  )


#### GO Plotting --------------

lung_func_enr_res_sig |> 
  separate(group,sep="-",
           into = c("exp_name","tf","direction"),
           remove = F) |> 
  group_by(group) |> 
  arrange(desc(-log10(p_adjust))) |> 
  slice_head(n=10) |> 
  ungroup() |> 
  mutate(term_name = tidytext::reorder_within(
    term_name,
    by = -log10(p_adjust),
    within = group)) |> 
  ggplot(aes(
    x=-log10(p_adjust+1e-30),
    y=term_name,
    fill=-log10(p_adjust+1e-30)
  ))+
  geom_col()+
  theme_custom()+
  tidytext::scale_y_reordered()+
  facet_wrap(~group,nrow=1,scales="free_y")

## Isoform-RNA Overlap ----------
library(ggvenn)

list("CD4 T cell RNA" = gsea_res_sig |> 
       filter(group == "CD4 T cell RNA") |> 
       pull(term_name),
     "CD4 T cell Isoforms" = enr_res_sig |> 
       filter(grepl("CD4 T cell Isoforms",group)) |> 
       pull(term_name)
) |> 
  ggvenn::ggvenn(
    fill_color = c("#A090A0FF","#F0E0F0FF"),
    stroke_size = .1
  )


list("CD16 Monocyte RNA" = gsea_res_sig |> 
       filter(group == "CD16 Monocyte RNA") |> 
       pull(term_name),
     "CD16 Monocyte Isoforms" = enr_res_sig |> 
       filter(grepl("CD16 Monocyte Isoforms",group)) |> 
       pull(term_name)
) |> 
  ggvenn::ggvenn(
    fill_color = c("#586888FF","#C8E0F8FF"),
    stroke_size = .1)

## Th-GSEA association ---------
se_long <- decon |> 
  tidybulk::scale_abundance(method = "TMM") |> 
  pivot_se() |> 
  mutate(log=log2(counts_scaled + 1))


paths <- gsea_res_sig |> 
  filter(group =="CD4 T cell RNA") |> 
  separate_rows(leadingEdge,sep = ",") |> 
  (\(df) split(df$leadingEdge,df$term_name))()


path_res <- map(names(paths),~{
  df <- se_long |> 
    filter(gene_name %in% paths[[.x]]) |> 
    group_by(.sample) |> 
    reframe(
      med_log=median(log),
      path=.x
    ) |> 
    inner_join(model_df,
               by=".sample")
  
  
  map_dfr(
    colnames(clr_mat),
    function(clr_col) {
      lm(
        as.formula(paste(
          clr_col,
          "~",
          "med_log + gli_age + gli_sex + fis"
        )),
        data = df
      ) |>
        tidy() |>
        filter(term == "med_log") |>
        mutate(cell = sub("^clr_", "", clr_col)) |>
        mutate(path =.x)
    }
  )
  
}) |> 
  bind_rows() |> 
  inner_join(
    gsea_res_sig |> 
      filter(term_name %in% names(paths)) |> 
      mutate(gsea_direction=ifelse(NES>0,"Up","Down")) |> 
      dplyr::select(term_name, gsea_direction),
    by = c("path" ="term_name")
  ) |> 
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

### Path Res Plotting -----------
path_res |> 
  group_by(cell_clean,gsea_direction) |> 
  reframe(
    med_est = median(estimate)
  ) |> 
  ggplot(aes(
    x=reorder(gsea_direction,med_est),
    y=reorder(cell_clean,med_est),
    fill=med_est
  ))+
  geom_tile()+
  scale_fill_gradient2(
    low="#5399b0",
    mid="white",
    high="#9B6981FF",
    midpoint = 0
  )+
  theme_custom()+
  labs(
    y="",
    x="GSEA Direction",
    fill="Median Estimate"
  )

## Save Data -----------
saveRDS(go_term_data,file="./results/useful_data/go_term_data.rds")

saveRDS(gsea_res_sig,file="./results/da_res/gsea_res_sig.rds")
saveRDS(enr_res_sig,file="./results/da_res/enr_res_sig.rds")
saveRDS(decon,file = "./results/da_res/decon.rds")

saveRDS(da_res_total,file = "./results/da_res/da_res_total.rds")

saveRDS(lung_func_enr_res_sig,file = "./results/da_res/lung_func_enr_res_sig.rds")
