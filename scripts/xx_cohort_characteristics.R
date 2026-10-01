## Categorical Plots -----------------


{
  p1 <- fev1_fvc_expom |> 
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
  
  
  p2 <- fev1_fvc_expom |> 
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
  
  
  p3 <- fev1_fvc_expom |> 
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
  
  p1 <- fev1_fvc_expom |> 
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

p2 <- fev1_fvc_expom |> 
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


p3 <- fev1_fvc_expom |> 
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

## PCA ---------------------

# MultiAssayExperiment Object
fev1_fvc_expom <- readRDS("./results/exwas/fev1_fvc_expom.rds")

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



rna_assays <- c("CD4 T cell RNA", "CD16 Monocyte RNA",
                "CD4 T cell Isoforms", "CD16 Monocyte Isoforms")

fev1_fvc_expom <- rna_assays |>
  purrr::reduce(normalize_rna_assay, .init = fev1_fvc_expom)

### plot pca --------------
fev1_fvc_expom@metadata$quality_control$pca = NULL

pca_res <- fev1_fvc_expom |>
  run_pca(log_trans_omics = TRUE,
          action = "get")



fev1_fvc_expom |>
  tidyexposomics::pivot_sample() |>
  ggplot(aes(x=PC1,
             y=PC2
             #color=pftfev1fvc_actual
             ))+
  geom_point(color="#1d2b64")+
  theme_custom()+
  # scale_color_gradientn(colors=c(
  #   "#000004FF","#B63679FF","#FCFDBFFF"
  # ))+
  labs(
    title="Sample PCA"
    #color=expression("FEV"[1]*"/FVC")
  )
