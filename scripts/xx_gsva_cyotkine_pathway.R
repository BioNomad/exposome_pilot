## General Pathways ---------------
source("./scripts/run_gsva_profiles.R")
source("./scripts/get_gsva_scores.R")
source("./scripts/pivot_se.R")

igs_names <- c(
  "GOBP_LEUKOCYTE_CHEMOTAXIS",
  "GOBP_CYTOKINE_PRODUCTION",
  "GOBP_T_CELL_TOLERANCE_INDUCTION",
  "GOBP_CELLULAR_SENESCENCE",
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_CELL_MATRIX_ADHESION",
  "GOBP_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL",
  "GOBP_T_CELL_DIFFERENTIATION",
  "GOBP_T_CELL_CHEMOTAXIS",
  "GOBP_T_CELL_CYTOKINE_PRODUCTION",
  "GOBP_MONOCYTE_CHEMOTAXIS",
  "GOBP_MONOCYTE_ACTIVATION",
  "GOBP_MONOCYTE_EXTRAVASATION",
  "GOBP_MONOCYTE_DIFFERENTIATION"
  
)

cd4_igs <- run_gsva_profiles(
  gs_names = igs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD4 T cell RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

cd16_igs <- run_gsva_profiles(
  gs_names = igs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD16 Monocyte RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)


## Cytokine GSVA -------------------
cytokine_gs_names <- msig |> 
  filter(grepl("_PRODUCTION",gs_name)) |> 
  filter(!grepl("REGULATION",gs_name)) |> 
  filter(grepl("^GOBP",gs_name)) |> 
  dplyr::count(gs_name) |> 
  pull(gs_name)


cd4_cytokines <- run_gsva_profiles(
  gs_names = cytokine_gs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD4 T cell RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

cd16_cytokines <- run_gsva_profiles(
  gs_names = cytokine_gs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD16 Monocyte RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

## Adhesion GSVA -------------------
adhesion_gs_names <- c(
  # Cell-matrix / ECM
  "GOBP_CELL_MATRIX_ADHESION",
  "GOBP_CELL_SUBSTRATE_ADHESION",
  "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_CELL_ADHESION_MEDIATED_BY_INTEGRIN",
  "GOBP_FOCAL_ADHESION_ASSEMBLY",
  "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT",
  "GOMF_INTEGRIN_BINDING",
  # Cell-cell / junctions
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_HOMOPHILIC_CELL_CELL_ADHESION",
  "GOBP_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN",
  "GOBP_TIGHT_JUNCTION_ORGANIZATION",
  "GOBP_EPITHELIAL_CELL_CELL_ADHESION",
  "GOMF_CADHERIN_BINDING",
  # Leukocyte / inflammatory
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL"
)

adhesion_gs_names <- msig |> filter(grepl("ADHESION|INTEGRIN|CADHERIN|SELECTIN|ADHERENS|TIGHT_JUNCTION|EXTRACELLULAR",gs_name)) |>
  dplyr::count(gs_name) |> 
  filter(n>10) |> 
  pull(gs_name)


cd4_ad <- run_gsva_profiles(
  gs_names = adhesion_gs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD4 T cell RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

cd16_ad <- run_gsva_profiles(
  gs_names = adhesion_gs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD16 Monocyte RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

## Mitochondria GSVA -------------------
mito_gs_names <- msig |> filter(grepl("MITOCHONDRIA",gs_name)) |>
  dplyr::count(gs_name) |> 
  filter(n>10) |> 
  pull(gs_name)


cd4_mito <- run_gsva_profiles(
  gs_names = mito_gs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD4 T cell RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

cd16_mito <- run_gsva_profiles(
  gs_names = mito_gs_names,
  expom = fev1_fvc_expom,
  exp_name = "CD16 Monocyte RNA",
  outcome = "pftfev1fvc_actual",
  covariates = c("gli_age" , "gli_sex" , "fis")
)

### Plotting General Pathways ---------------

cd4_igs$lm_res |> 
  bind_rows(
    cd16_igs$lm_res 
  ) |> 
  ggplot(aes(
    x = reorder(exp_name, estimate),
    y = reorder(term, estimate),
    fill = estimate
  )) +
  geom_tile() +
  geom_point(
    data = ~ filter(.x, p.value < 0.1),
    aes(shape = "P < 0.1"),
    fill = "grey90",
    size = 3
  ) +
  theme_custom() +
  scale_fill_gradient2(
    low = "#5399b0",
    mid = "white",
    high = "#9B6981FF",
    midpoint = 0
  ) +
  scale_shape_manual(
    values = c("P < 0.1" = 22),
    name = ""
  ) +
  scale_x_discrete(labels = c(
    "CD4 T cell RNA"    = expression(CD4^"+" ~ "T cell"),
    "CD16 Monocyte RNA" = expression(CD16^"+" ~ "Monocyte")
  )) +
  labs(
    x = "",
    y = "",
    fill = "Estimate"
  ) +
  theme(
    axis.text.x = element_text(angle = 90,vjust =.5),
    axis.text.y = element_text(size = 11)
  )
### Plotting Cytokines ---------------

cd4_cytokines$lm_res |> 
  bind_rows(
    cd16_cytokines$lm_res 
  ) |> 
  filter(!is.na(p.value)) |> 
  ggplot(aes(
    x = reorder(exp_name, estimate),
    y = reorder(term, estimate),
    fill = estimate
  )) +
  geom_tile() +
  geom_point(
    data = ~ filter(.x, p.value < 0.1),
    aes(shape = "P < 0.1"),
    fill = "grey90",
    size = 3
  ) +
  theme_custom() +
  scale_fill_gradient2(
    low = "#5399b0",
    mid = "white",
    high = "#9B6981FF",
    midpoint = 0
  ) +
  scale_shape_manual(
    values = c("P < 0.1" = 22),
    name = ""
  ) +
  scale_x_discrete(labels = c(
    "CD4 T cell RNA"    = expression(CD4^"+" ~ "T cell"),
    "CD16 Monocyte RNA" = expression(CD16^"+" ~ "Monocyte")
  )) +
  labs(
    x = "",
    y = "",
    fill = "Estimate"
  ) +
  theme(
    axis.text.x = element_text(angle = 90,vjust =.5),
    axis.text.y = element_text(size = 11)
  )

### Plotting Adhesion ------------------
cd4_ad$lm_res |> 
  bind_rows(
    cd16_ad$lm_res 
  ) |> 
  ggplot(aes(
    x = reorder(exp_name, estimate),
    y = reorder(term, estimate),
    fill = estimate
  )) +
  geom_tile() +
  geom_point(
    data = ~ filter(.x, p.value < 0.1),
    aes(shape = "P < 0.1"),
    fill = "grey90",
    size = 3
  ) +
  theme_custom() +
  scale_fill_gradient2(
    low = "#5399b0",
    mid = "white",
    high = "#9B6981FF",
    midpoint = 0
  ) +
  scale_shape_manual(
    values = c("P < 0.1" = 22),
    name = ""
  ) +
  scale_x_discrete(labels = c(
    "CD4 T cell RNA"    = expression(CD4^"+" ~ "T cell"),
    "CD16 Monocyte RNA" = expression(CD16^"+" ~ "Monocyte")
  )) +
  labs(
    x = "",
    y = "",
    fill = "Estimate"
  ) +
  theme(
    axis.text.x = element_text(angle = 90,vjust =.5),
    axis.text.y = element_text(size = 11)
  )
### Focused Cytokine Results -------------
cd4_cytokines$lm_res |> 
  filter(p.value<0.1) |> 
  filter(term != "GOBP_MYELOID_DENDRITIC_CELL_CYTOKINE_PRODUCTION") |> 
  bind_rows(
    cd16_cytokines$lm_res |> 
      filter(p.value<0.1) |> 
      filter(term != "GOBP_CYTOKINE_PRODUCTION")
  ) |> 
  mutate(term = gsub("GOBP_|_PRODUCTION|","",term)) |> 
  # mutate(term = gsub("INTERLEUKIN_","IL",term)) |> 
  # mutate(term = case_when(
  #   term == "MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA" ~ "MIP1A",
  #   term == "CHEMOKINE_C_C_MOTIF_LIGAND_5" ~ "CCL5",
  #   term == "NEUROTROPHIN" ~ "Neurotrophin",
  #   term == "TYPE_II_INTERFERON" ~ "Type II IFN",
  #   term == "TYPE_III_INTERFERON" ~ "Type III IFN",
  #   .default = term
  # )) |> 
  ggplot(aes(
    x=abs(estimate),
    y=reorder(term,abs(estimate)),
    fill=estimate
  ))+
  geom_col()+
  theme_custom()+
  facet_grid(
    exp_name ~ .,
    space = "free",
    scales = "free",
    labeller = as_labeller(c(
      "CD4 T cell RNA"    = "CD4\u207a T cell",
      "CD16 Monocyte RNA" = "CD16\u207a Monocyte"
    ))
  )+
  scale_fill_gradient2(
    low = "#5399b0",
    mid = "white",
    high = "#9B6981FF",
    midpoint = 0
  ) +
  labs(
    x="|Estimate|",
    y="",
    fill="Estimate"
  )+
  scale_x_continuous(
    breaks = c(0,.2,.4)
  )+
  theme(
    axis.text.x = element_text(angle=0,hjust=.3),
    #axis.text.y = element_text(face="italic")
  )

## Adhesion Distributions --------------
cd4_ad$lm_res |> 
  bind_rows(
    cd16_ad$lm_res 
  ) |> 
  mutate(exp_name=factor(exp_name,levels=c(
    "CD4 T cell RNA",
    "CD16 Monocyte RNA"
  ))) |> 
  ggplot(aes(
    x=exp_name,
    y=estimate,
    fill=exp_name
  ))+
  geom_boxplot(alpha=.5)+
  theme_custom()+
  geom_hline(yintercept = 0,linetype="dashed")+
  stat_compare_means(
    label = "p.signif",ref.group = "CD4 T cell RNA"
  )+
  guides(fill="none")+
  scale_y_continuous(expand = c(0,.1))+
  labs(
    x="",
    y="Adhesion Estimate"
  )+
  scale_fill_manual(values=c(
    "CD16 Monocyte RNA"="#C2697FFF",
    "CD4 T cell RNA"="#B4B9E0FF"
  ))+
  scale_x_discrete(labels = c(
    "CD4 T cell RNA"    = expression(CD4^"+" ~ "T cell"),
    "CD16 Monocyte RNA" = expression(CD16^"+" ~ "Monocyte")
  )) 

## Immunoglobulins -----------------
da_res_total |>
  filter(exp_name == "Proteomics") |> 
  mutate(ig = case_when(
    grepl("^IG",feature_id) ~ "Immunoglobulins",
    .default = "All Other Proteins"
  )) |> 
  ggplot(aes(
    x=ig,
    y=logfc,
    fill=ig
  ))+
  geom_boxplot(alpha=.5)+
  theme_custom()+
  geom_hline(yintercept = 0,linetype="dashed")+
  stat_compare_means(
    label = "p.signif",
    ref.group = "All Other Proteins"
  )+
  labs(
    x="",
    y=expression("Log"[2]*"FC")
  )+
  scale_y_continuous(expand = c(0,1))+
  guides(fill="none")+
  scale_fill_manual(values=c(
    "All Other Proteins"="grey90",
    "Immunoglobulins"="#3E6992FF"
  ))

## Th Correlations with Cytokines ----------------
cytokine_terms <- cd4_cytokines$lm_res |> 
  filter(p.value < 0.1) |> 
  pull(term)

cd4_cytokines$profiles |> 
  dplyr::select(.sample, 
                all_of(cytokine_terms)) |> 
  inner_join(
    tidybulk::pivot_sample(decon) |> 
      dplyr::select(.sample, 
                    all_of(cell_cols)),
    by = ".sample"
  ) |> 
  pivot_longer(
    all_of(cytokine_terms),
    names_to = "cytokine",
    values_to = "cytokine_score") |>
  pivot_longer(all_of(cell_cols), 
               names_to = "cell_type",
               values_to = "proportion") |>
  summarise(
    estimate = cor(cytokine_score, 
                   proportion,
                   method = "spearman"),
    p.value  = cor.test(cytokine_score,
                        proportion,
                        method = "spearman")$p.value,
    .by = c(cytokine, cell_type)
  ) |>
  mutate(cytokine = gsub("GOBP_|_PRODUCTION|","",cytokine)) |> 
  # mutate(cytokine = gsub("INTERLEUKIN_","IL",cytokine)) |> 
  # mutate(cytokine = case_when(
  #   cytokine == "MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA" ~ "MIP1A",
  #   cytokine == "CHEMOKINE_C_C_MOTIF_LIGAND_5" ~ "CCL5",
  #   cytokine == "NEUROTROPHIN" ~ "Neurotrophin",
  #   cytokine == "TYPE_II_INTERFERON" ~ "Type II IFN",
  #   cytokine == "TYPE_III_INTERFERON" ~ "Type III IFN",
  #   .default = cytokine
  # )) |> 
  ggplot(aes(
    x = reorder(cell_type,estimate),
    y = reorder(cytokine,estimate),
    fill = estimate
  )) +
  geom_tile() +
  geom_point(
    data = ~ filter(.x, p.value < 0.1),
    aes(shape = "P < 0.1"),
    fill  = "grey90",
    size  = 3
  ) +
  scale_shape_manual(
    values = c("P < 0.1" = 22),
    name = "") +
  scale_fill_gradient2(
    high  = "#C33764",
    mid   = "white",
    low   = "#1D2671",
    midpoint = 0
  ) +
  scale_x_discrete(
    labels = \(x) str_remove(x, "Memory\\.CD4\\.T\\.cell\\.") |>
      str_replace_all("\\.", " ")) +
  theme_custom() +
  theme(
    # axis.text.x = element_text(angle = 90, vjust = .5)
    axis.text.x = element_text(angle = 60, just = 1)
  ) +
  labs(x = "", 
       y = "",
       fill = expression(~rho))
