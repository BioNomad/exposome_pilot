## Load Libraries ------------

library(tidyverse)
library(janitor)
library(ggpubr)
library(tidyexposomics)
library(MultiAssayExperiment)

source("./scripts/bin/internals.R")
source("./scripts/bin/labels_clean.R")
## Load Data --------------
# Meta Data
aw_meta_cv2 <- readRDS("./results/input_data/aw_meta_cv2.rds")

# Codebook
aw_cb <- readRDS("./results/input_data/aw_cb.rds")

# Abundance/Expression Data
cd4_gene_counts   <- readRDS("./results/input_data/cd4_gene_counts.rds")
cd16_gene_counts  <- readRDS("./results/input_data/cd16_gene_counts.rds")

cd4_isoform_counts  <- readRDS("./results/input_data/cd4_isoform_counts.rds")
cd16_isoform_counts <- readRDS("./results/input_data/cd16_isoform_counts.rds")

prot_abd <- readRDS("./results/input_data/prot_abd.rds")
adduct   <- readRDS("./results/input_data/adduct.rds")
adduct_load   <- readRDS("./results/input_data/adduct_load.rds")

# Feature Data
gene_fdata    <- readRDS("./results/input_data/gene_fdata.rds")
isoform_fdata <- readRDS("./results/input_data/isoform_fdata.rds")
prot_fdata    <- readRDS("./results/input_data/prot_fdata.rds")
adduct_fdata  <- readRDS("./results/input_data/adduct_fdata.rds")
adduct_load_fdata <- readRDS("./results/input_data/adduct_load_fdata.rds")

## Create the omics list -----------
omics_list <- list(
  "CD4 T cell RNA"       = cd4_gene_counts,
  "CD16 Monocyte RNA"      = cd16_gene_counts,
  "CD4 T cell Isoforms"  = cd4_isoform_counts,
  "CD16 Monocyte Isoforms" = cd16_isoform_counts,
  "Proteomics"            = prot_abd,
  "Adductomics"           = adduct,
  "Adduct Burden"      = adduct_load
)

fdata <- list(
  "CD4 T cell RNA"       = gene_fdata,
  "CD16 Monocyte RNA"      = gene_fdata,     
  "CD4 T cell Isoforms"  = isoform_fdata,
  "CD16 Monocyte Isoforms" = isoform_fdata,   
  "Proteomics"            = prot_fdata,
  "Adductomics"           = adduct_fdata,
  "Adduct Burden"      = adduct_load_fdata
)

rm(list = ls()[grepl("_fdata|counts|_abd|adduct", ls())])


expom <- create_exposomicset(
  codebook = aw_cb,
  exposure = aw_meta_cv2,
  omics = omics_list,
  row_data = fdata
)

# Keep samples with sex and age
expom <- expom[,!is.na(expom$age) & !is.na(expom$gli_sex)]

exp_cat <- c(
  "Indoor Air",
  "Parabens",
  "Phenols",
  "Antimicrobials",
  
  "Dust Allergens",
  "IgE Respiratory",
  "IgE Superantigen",
  #"IgE Total",
  "IgE Food",
  
  "Serum Essential Metals",
  "Serum Non-Essential Metals",
  "Urine Essential Metals",
  "Urine Non-Essential Metals"
  
  )

exp_vars <- extract_results(exposomicset = expom,result = "codebook") |>
  filter(category %in% exp_cat) |> 
  pull(variable)

# Setting colors for exposure categories
cat_colors <- scales::alpha(get_palette("aaas",13), 0.5) |> 
  (\(colors){
    names(colors) <- aw_cb |> 
      filter(variable %in% exp_vars) |> 
      pull(category) |> 
      unique()
    colors
  })()

## ----------------
expom |>
    plot_missing(
        plot_type = "summary",
        threshold = 20
    )


## Exposure LOD Imputation -------------
# Grab codebook variables that are actually in colData
codebook_vars <- extract_results(expom, result = "codebook") |>
  filter(variable %in% (
    expom |> tidyexposomics::pivot_sample() |> colnames()
  ))

# LOD-style imputation
lod_cats <- c(
  
  "Indoor Air",
  "Parabens",
  "Phenols",
  "Antimicrobials",
  
  "Dust Allergens",
  "IgE Respiratory",
  "IgE Superantigen",
  "IgE Total",
  "IgE Food",
  
  "Serum Essential Metals",
  "Serum Non-Essential Metals",
  "Urine Essential Metals",
  "Urine Non-Essential Metals",
  
  "Serum Cytokines",
  "Nasal Cytokines"
  
)

# Split into two lists of variables
lod_vars <- codebook_vars |> 
  filter(category %in% lod_cats) |> 
  pull(variable)

# remove variables with no variance
lod_vars <- expom |>
  tidyexposomics::pivot_sample() |>
  dplyr::select(all_of(lod_vars)) |>
  apply(2,function(col){
    var <- col[!is.na(col)] |> var();
    return(var>0)
  }) |> 
  (\(bool) bool[bool %in% c(TRUE)])() |> 
  names()

# now confirm the missingness of each:
lod_vars <- expom |> 
  tidyexposomics::pivot_sample() |> 
  dplyr::select(all_of(lod_vars)) |> 
  naniar::miss_var_summary() |> 
  filter(pct_miss!=100) |> 
  pull(variable)

# looks like half the samples were just not profiled
# for some chemicals, removing these since 
# they are not non-detects
lod_vars <- lod_vars[!grepl("_sg",lod_vars)]
  
# now let's impute
expom <- expom |>
  run_impute_missing(
    exposure_impute_method = "lod_sqrt2",
    exposure_cols = lod_vars
  ) 


## Filter Missing Exposures/Omics Features -----------


expom <- expom |> 
  filter_missing(na_thresh = 20)

# update exp vars after filtering
exp_vars <- exp_vars[
  exp_vars %in% names(colData(expom))
]



## missForest Imputation --------------
# Define missForest variables to impute
mf_vars <- expom |>
    tidyexposomics::pivot_sample() |> 
  naniar::miss_var_summary() |> 
  filter(pct_miss>0) |> 
  pull(variable)

# Now run imputation by bucket
expom <- expom |>
  run_impute_missing(
    exposure_impute_method = "missforest",
    exposure_cols = mf_vars,
    omics_impute_method  = "missforest",
    omics_to_impute = c(
      "Proteomics","Adductomics")
  ) 
  


## Filtering Omics --------------------------
expom <- expom |>
    filter_omics(
        method = "expression",
        assays = "CD4 T cell RNA",
        assay_name = 1,
        min_value = 1,
        min_prop = 0.3
    ) |>
  filter_omics(
        method = "expression",
        assays = "CD16 Monocyte RNA",
        assay_name = 1,
        min_value = 1,
        min_prop = 0.3
    ) |>
  filter_omics(
        method = "expression",
        assays = "CD4 T cell Isoforms",
        assay_name = 1,
        min_value = 1,
        min_prop = 0.3
    ) |>
  filter_omics(
        method = "expression",
        assays = "CD16 Monocyte Isoforms",
        assay_name = 1,
        min_value = 1,
        min_prop = 0.3
    ) |>
    filter_omics(
        method = "variance",
        assays = "Proteomics",
        assay_name = 1,
        min_var = 0.01
    ) 


## Normality/Transform Exposures ----------------------------

num_exp_vars <- expom |>
  tidyexposomics::pivot_sample() |> 
  dplyr::select(all_of(exp_vars)) |> 
  dplyr::select_if(is.numeric) |> 
  dplyr::select_if(~length(unique(.))>3) |> 
  # keeping temperature as is
  dplyr::select(-c(iv_tempavg)) |> 
  colnames()

# Check variable normality & transform variables
expom <- expom |>
    # Check variable normality
    run_normality_check(action = "add")  |> 
    # Transform variables
    transform_exposure(
        transform_method = "log2",
        exposure_cols = num_exp_vars
    )


## PCA -----------------

# Perform principal component analysis
expom <- expom |>
    run_pca(
        log_trans_exp = TRUE,
        log_trans_omics = TRUE,
        action = "add"
    )
# Plot principal component analysis results
# expom |>
#     plot_pca()

## Removing Outliers ----------------------------
#Filter out sample outliers
# expom <- expom |>
#   filter_sample_outliers(
#     outliers = c("s937", "s3331","s3692")
#   )


## High Exposure Analysis ------------

meta <- expom |> 
  pivot_sample()

quants <- map(num_exp_vars,~{
  quants <- quantile(meta[[.x]])
  data.frame(
    variable=.x,
    quant_25 = quants["25%"],
    quant_50 = quants["50%"],
    quant_75 = quants["75%"]
  )
}) |> 
  bind_rows() 

expom |> 
  pivot_sample() |> 
  pivot_longer(
    cols = num_exp_vars,
    names_to = "variable",
    values_to = "level"
  ) |> 
  inner_join(
    quants ,
    by ="variable"
  ) |> 
  mutate(over=level>quant_75) |> 
  dplyr::select(variable,.sample,level,over,quant_75) |> 
  filter(over==TRUE) |> 
  dplyr::count(.sample) |> 
  view()

expom |> 
  pivot_sample() |> 
  pivot_longer(
    cols = num_exp_vars,
    names_to = "variable",
    values_to = "level"
  ) |> 
  inner_join(
    quants ,
    by ="variable"
  ) |> 
  mutate(over=level>quant_75) |> 
  dplyr::select(variable,.sample,level,over,quant_75) |> 
  filter(over==TRUE) |> 
  dplyr::count(.sample) |> 
  ggplot(aes(x=n))+
  geom_density(fill="#474554")+
  theme_custom()+
  labs(
    x=expression("No. of Exposures"[High]),
    y="Density"
  )

## Association with Asthma ------------

expom |> 
  pivot_sample() |> 
  pivot_longer(
    cols = num_exp_vars,
    names_to = "variable",
    values_to = "level"
  ) |> 
  inner_join(
    quants ,
    by ="variable"
  ) |> 
  mutate(over=level>quant_75) |> 
  dplyr::select(variable,.sample,level,over,quant_75) |> 
  filter(over==TRUE) |> 
  dplyr::count(.sample) |> 
  inner_join(
    expom |> 
      pivot_sample(),
    by=".sample"
  ) |> 
  ggplot(aes(
    x=pftfev1fvc_actual,
    y=n
  ))+
  geom_point(color="black",size=2.5,shape=21,fill="grey50")+
  theme_custom()+
  stat_cor(p.accuracy  = 0.001,
           method = "spearman",
           label.y =32)+
  stat_smooth(method=lm,color="#AD5A6BFF")+
  labs(
    x=expression("FEV"[1]*"/FVC"),
    y=expression("No. of Exposures"[High])
  )


## Exposure Variablity ---------------


expom |> 
  run_summarize_exposures(
    exposure_cols = num_exp_vars,
    action = "get") |> 
  arrange(desc(coef_var)) |> 
  #slice_head(n=10) |> 
  mutate() |> 
  ggplot(aes(
    x=coef_var,
    y=reorder(clean_name,coef_var)
  ))+
  geom_col()+
  geom_vline(xintercept = .5,linetype="dashed")+
  theme_custom()


expom |> 
  run_summarize_exposures(
    exposure_cols = num_exp_vars,
    action = "get") |> 
  ggplot(aes(
    x=reorder(category,coef_var),
    y=coef_var,
    fill=category,
    color=category
  ))+
  geom_jitter(alpha=.1)+
  geom_boxplot(color="black")+
  theme_custom()+
  theme(axis.text.x=element_text(angle=70))+
  geom_hline(yintercept = 0.5,
             linetype="dashed")+
  guides(color="none",
         fill="none")+
  scale_color_manual(values=cat_colors)+
  scale_fill_manual(values=cat_colors)+
  labs(
    x="",
    y="Coef. Var."
  )+
  coord_flip()+
  theme(axis.text.x=element_text(angle = 0,hjust=.3))+
  scale_y_continuous(expand=c(.2,0))

## Correlation of Exposures with joint omics space--------------

expom <- expom |>
    run_correlation(
        feature_type = "pcs",
        exposure_cols = num_exp_vars,
        n_pcs = 30,
        action = "add",
        correlation_cutoff = 0,
        pval_cutoff = 1
    )

expom |>
    plot_correlation_tile(
        feature_type = "pcs",
        pval_cutoff = 0.1
    )

## Categorical Plots -----------------


{
  p1 <- expom |> 
    tidyexposomics::pivot_sample() |> 
    mutate(sex_lab=case_when(
      gli_sex == "M" ~ "Male",
      gli_sex == "F" ~ "Female"
    )) |> 
    make_donut(fill_var = "sex_lab",
               fill_vals = c(
                 "Male"="#365C83FF",
                 "Female"="#F6B6C2FF"
               ),
               title="Sex")
  
  
  p2 <- expom |> 
    tidyexposomics::pivot_sample() |> 
    mutate(anc_lab=case_when(
      black == "black" ~ "Black",
      black == "non-black" ~ "Other Ancestry"
    )) |> 
    make_donut(fill_var = "anc_lab",
               fill_vals = c(
                 "Black"="#A89F8EFF",
                 "Other Ancestry"="#7B7987FF"
               ),
               title="Ancestry")
  
  
  p3 <- expom |> 
    tidyexposomics::pivot_sample() |> 
    mutate(income=case_when(
      income5 == "RA/DK" ~ "Not Reported",
      .default = income5
    )) |> 
    make_donut(fill_var = "income",
               fill_vals = c(
                 "<$15k" = "#E3CACFFF",
                 "$15k-$35k" = "#cfb8d6",
                 "$35k-50k" = "#af9cb5",
                 ">=50k" = "#453d47",
                 "Not Reported"="grey90"
               ),
               title="Household Income")
  
  # (p1 | p2 | p3)
  (p1 / p2 / p3)
}


# Continuous Plots ---------------
{
  
  p1 <- expom |> 
    tidyexposomics::pivot_sample() |> 
    ggplot(aes(
      x=gli_age
    ))+
    geom_density(fill="#1D2731FF" |> scales::alpha(.8))+
    theme_custom()+
    labs(
      x="Age",
      #title = "Age",
      y = "Density"
    )+
    theme(
      axis.text.x = element_text(angle=0,hjust=0.5)
    )
  
  p2 <- expom |> 
    tidyexposomics::pivot_sample() |> 
    ggplot(aes(
      x=bmizscore
    ))+
    geom_density(fill="#1D2731FF" |> scales::alpha(.8))+
    theme_custom()+
    labs(
      x="Z-scored BMI",
      #title = "Z-scored BMI",
      y = "Density"
    )+
    theme(
      axis.text.x = element_text(angle=0,hjust=0.5)
    )
  
  
  p3 <- expom |> 
    tidyexposomics::pivot_sample() |> 
    ggplot(aes(
      x=pftfev1fvc_actual
    ))+
    geom_density(fill="#1D2731FF" |> scales::alpha(.8))+
    theme_custom()+
    labs(
      x=expression("FEV"[1]*"/FVC"),
      #title = bquote(bold("FEV")[1] ~ bold("/FVC")),
      y = "Density"
    )+
    theme(
      axis.text.x = element_text(angle=0,hjust=0.5)
    )
  
  #(p1 | p2 | p3)
  (p1 / p2 / p3)
}

## Example LM Code -------------
iris |> 
  filter(Petal.Length>3) |> 
  ggplot(aes(
    x=Petal.Length,
    y=Petal.Width
  ))+
 geom_point(color="black",size=2,shape=21,fill="grey50")+
  stat_smooth(method=lm,color="#AD5A6BFF")+
  theme_custom()+
  theme(
    axis.text = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(
    x="Median Ledge\nGene Expression",
    y="Th Cell"
  )
## Save Data ------------------------------
saveRDS(expom,file="./results/qc/expom.rds")
saveRDS(exp_vars,file="./results/qc/exp_vars.rds")
saveRDS(num_exp_vars,file="./results/qc/num_exp_vars.rds")

