# --- Run Limma -------------
run_limma <- function(
    se,
    formula,
    levels,
    abundance_col="counts",
    study=NULL,
    subset_filters = list(),
    method = c("trend", "voom"),
    contrasts = NULL,
    robust = TRUE,
    log_trans = F,
    minimum_counts = NULL,
    minimum_proportion = NULL
    
) {
  message("Study: ",study)
  method <- match.arg(method)
  
  # Subset to study
  if(!is.null(study)){
    se <- se[, grepl(study, se$study)]
  }
  
  # Apply optional filters
  for (flt in subset_filters) {
    se <- se[, grepl(flt$pattern, 
                     colData(se)[[flt$col]])]
  }
  
  # Set levels
  if(!is.null(levels)){
    se[[(all.vars(formula)[1])]] <- factor(
      se[[(all.vars(formula)[1])]],
      levels = levels
    )
  }
  message("Tissue: ",unique(se$tissue))
  
  
  # Filter the SE
  if(!is.null(minimum_counts) & !is.null(minimum_proportion)){
    se <- se |>
      identify_abundant(
        minimum_counts = minimum_counts,
        minimum_proportion = minimum_proportion
      ) |> 
      keep_abundant() 
  }
  
  # Extract expression matrix
  mat <- SummarizedExperiment::assay(se, abundance_col)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  
  # Remove NA values
  if(any(is.na(mat))) mat <- na.omit(mat)
  
  message(paste(dim(mat),collapse = ","))
  
  # If negative add the pseudo-count
  if(min(mat)<0){
    mat <- mat + abs(min(mat,na.rm = T))
  }
  
  # Log2-transform if needed
  if (log_trans) {
    mat <- log2(mat + 1)
  }
  
  # Extract metadata and create design matrix
  df <- as.data.frame(SummarizedExperiment::colData(se))
  message(unique(df$diagnosis))
  design <- stats::model.matrix(formula, data = df)
  
  message(unique(df$diagnosis))
  
  # Determine limma fitting method
  if (method == "voom") {
    message("Running limma-voom...")
    lib_sizes <- edgeR::calcNormFactors(edgeR::DGEList(counts = mat))$samples$lib.size
    v <- limma::voom(mat, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
  } else {
    message("Running limma-trend...")
    fit <- limma::lmFit(mat, design)
  }
  
  fit <- limma::eBayes(fit, trend = (method == "trend"), robust = robust)
  
  # Compute contrasts or coefficients
  results <- list()
  if (!is.null(contrasts)) {
    cm <- limma::makeContrasts(contrasts = contrasts, levels = design)
    fit2 <- limma::contrasts.fit(fit, cm)
    fit2 <- limma::eBayes(fit2, trend = (method == "trend"), robust = robust)
    for (i in seq_len(ncol(cm))) {
      tb <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none")
      tb$contrast <- colnames(cm)[i]
      results[[i]] <- tb
    }
  } else {
    coef_idx <- setdiff(seq_len(ncol(design)), which(colnames(design) == "(Intercept)"))
    for (k in coef_idx) {
      tb <- limma::topTable(fit, coef = k, number = Inf, sort.by = "none")
      tb$contrast <- colnames(design)[k]
      results[[length(results) + 1]] <- tb
    }
  }
  
  out <- dplyr::bind_rows(results)
  
  # Add annotations and metadata
  out <- out |>
    tibble::rownames_to_column("feature") |>
    dplyr::mutate(
      feature = gsub("\\.\\.\\..*", "", feature),
      method = paste0("limma_", method)
    ) |>
    dplyr::left_join(
      se |> tidybulk::pivot_transcript(),
      by = c("feature" = ".feature")
    ) |>
    dplyr::filter(grepl(all.vars(formula)[1], contrast)) |>
    tibble::as_tibble()
  
  return(out)
}

# --- Suggest Covariates ---------
# Identify potential covariates by study/tissue via PC associations
# Usage:
#   suggest_covariates(se = gene_se, outcome = "diagnosis",
#                      group_vars = c("study","tissue"),
#                      abundance_col = "counts", log_trans = NULL)
# Apply to the full SE; grouping is handled internally so PCs are computed
# within each study × tissue subset.
# feasible:	basic structural validity	enough samples, non-constant	TRUE
# recommended: 	strong PC correlation	− log10(p) ≥ 2	TRUE
# vif_ok:	non-redundant	VIF < 5	TRUE
# outcome_p: 	not collinear with diagnosis	p > 0.05	TRUE
suggest_covariates <- function(
    se,
    outcome = "diagnosis",
    group_vars = c("study", "tissue"),
    abundance_col = NULL,
    log_trans = NULL,
    pc_k = 5,
    min_non_na_prop = 0.7,
    max_levels = 8,
    min_per_level = 3,
    top_var_features = 3000,
    candidate_cols = NULL,
    exclude_additional = NULL
) {
  df <- as.data.frame(SummarizedExperiment::colData(se))
  
  # Resolve grouping columns present
  group_vars <- intersect(group_vars, colnames(df))
  if (length(group_vars) == 0) {
    df$.__all__ <- "all"
    group_vars <- ".__all__"
  }
  
  # Resolve assay and transform defaults
  if (is.null(abundance_col)) abundance_col <- SummarizedExperiment::assayNames(se)[1]
  if (is.null(log_trans))      log_trans     <- identical(abundance_col, "counts")
  
  # Build groups (tidy)
  groups_df <- dplyr::distinct(df, dplyr::across(dplyr::all_of(group_vars)))
  
  # Utility: compute PCs for a subset of samples
  compute_pcs_for_subset <- function(se_sub) {
    mat <- SummarizedExperiment::assay(se_sub, abundance_col)
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    if (log_trans) mat <- log2(mat + 1)
    keep <- matrixStats::rowSds(mat, na.rm = TRUE) > 0
    mat <- mat[keep, , drop = FALSE]
    if (nrow(mat) == 0 || ncol(mat) < 3) return(NULL)
    if (nrow(mat) > top_var_features) {
      sds <- matrixStats::rowSds(mat)
      idx <- order(sds, decreasing = TRUE)[seq_len(top_var_features)]
      mat <- mat[idx, , drop = FALSE]
    }
    pc <- tryCatch(stats::prcomp(t(mat), center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (is.null(pc)) return(NULL)
    k <- min(pc_k, ncol(pc$x))
    pcs <- pc$x[, seq_len(k), drop = FALSE]
    colnames(pcs) <- paste0("PC", seq_len(k))
    pcs
  }
  
  # Exclusions for candidate covariates
  exclude_cols <- c(".sample","sra","gse_id","title","study","tissue","source", outcome)
  if (!is.null(exclude_additional)) exclude_cols <- unique(c(exclude_cols, exclude_additional))
  
  # Helper: association p-value with PCs
  assoc_min_p <- function(v, pcs_mat) {
    pvals <- apply(pcs_mat, 2, function(pc) {
      if (is.numeric(v)) {
        suppressWarnings(stats::cor.test(pc, v, method = "spearman")$p.value)
      } else {
        suppressWarnings(stats::anova(stats::lm(pc ~ as.factor(v)))$`Pr(>F)`[1])
      }
    })
    pvals <- as.numeric(pvals)
    pvals[!is.finite(pvals)] <- NA_real_
    suppressWarnings(min(pvals, na.rm = TRUE))
  }
  
  # Helper: safe outcome p-value (collinearity diagnostic)
  safe_p_out <- function(v, outcome_vec, outcome_fac) {
    outcome_has_var <- (is.numeric(outcome_vec) && is.finite(stats::sd(outcome_vec, na.rm = TRUE)) && stats::sd(outcome_vec, na.rm = TRUE) > 0) ||
      (!is.null(outcome_fac) && nlevels(outcome_fac) >= 2)
    if (!outcome_has_var) return(NA_real_)
    tryCatch({
      if (is.numeric(v) && !is.null(outcome_fac) && nlevels(outcome_fac) >= 2) {
        if (stats::sd(v, na.rm = TRUE) == 0) return(NA_real_)
        suppressWarnings(stats::anova(stats::lm(v ~ outcome_fac))$`Pr(>F)`[1])
      } else if (!is.numeric(v) && !is.null(outcome_fac) && nlevels(outcome_fac) >= 2) {
        fac_v <- droplevels(as.factor(v))
        tbl <- table(fac_v, outcome_fac, useNA = "no")
        if (all(dim(tbl) > 1)) suppressWarnings(stats::chisq.test(tbl)$p.value) else NA_real_
      } else if (is.numeric(v) && is.numeric(outcome_vec) &&
                 stats::sd(v, na.rm = TRUE) > 0 &&
                 stats::sd(outcome_vec, na.rm = TRUE) > 0) {
        suppressWarnings(stats::cor.test(v, outcome_vec, method = "spearman")$p.value)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
  }
  
  # Helper: compute VIFs for feasible covariates (per group)
  #   compute_vif_tbl <- function(df_sub, covar_cols, outcome_name) {
  #   if (!requireNamespace("car", quietly = TRUE)) {
  #     return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
  #   }
  #   cols <- unique(stats::na.omit(c(outcome_name, covar_cols)))
  #   dat <- df_sub[, cols, drop = FALSE]
  # 
  #   # Outcome to numeric
  #   outcome_num <- as.numeric(as.factor(dat[[outcome_name]]))
  #   dat$outcome_num <- outcome_num
  # 
  #   pred_cols <- setdiff(colnames(dat), c(outcome_name, "outcome_num"))
  #   if (length(pred_cols) < 2) {
  #     return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
  #   }
  # 
  #   dat_cc <- tidyr::drop_na(dat, dplyr::all_of(c("outcome_num", pred_cols)))
  #   if (nrow(dat_cc) < 4 || length(unique(dat_cc$outcome_num)) < 2) {
  #     return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
  #   }
  # 
  #   frm <- stats::as.formula(paste("outcome_num ~", paste(pred_cols, collapse = " + ")))
  #   fit <- tryCatch(stats::lm(frm, data = dat_cc), error = function(e) NULL)
  #   if (is.null(fit)) {
  #     return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
  #   }
  # 
  #   # Detect aliased coefficients and refit if needed
  #   aliased <- names(coef(fit))[is.na(coef(fit))]
  #   if (length(aliased) > 0) {
  #     valid_preds <- setdiff(pred_cols, aliased)
  #     if (length(valid_preds) < 2) {
  #       return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
  #     }
  #     frm2 <- stats::as.formula(paste("outcome_num ~", paste(valid_preds, collapse = " + ")))
  #     fit <- tryCatch(stats::lm(frm2, data = dat_cc), error = function(e) NULL)
  #   }
  # 
  #   v <- tryCatch(car::vif(fit), error = function(e) NULL)
  #   if (is.null(v)) {
  #     return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
  #   }
  # 
  #   # Handle both numeric and GVIF table outputs
  #   if (is.null(dim(v))) {
  #     tibble::tibble(variable = names(v), vif = as.numeric(v), gvif_adj = NA_real_)
  #   } else {
  #     v <- as.data.frame(v)
  #     v$variable <- rownames(v)
  #     v$gvif_adj <- v$GVIF^(1/(2 * v$Df))
  #     tibble::tibble(variable = v$variable, vif = NA_real_, gvif_adj = v$gvif_adj)
  #   }
  # }
  compute_vif_tbl <- function(df_sub, covar_cols, outcome_name) {
    if (!requireNamespace("car", quietly = TRUE)) {
      return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
    }
    
    cols <- unique(stats::na.omit(c(outcome_name, covar_cols)))
    dat <- df_sub[, cols, drop = FALSE]
    dat$outcome_num <- as.numeric(as.factor(dat[[outcome_name]]))
    
    pred_cols <- setdiff(colnames(dat), c(outcome_name, "outcome_num"))
    if (length(pred_cols) < 2) {
      return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
    }
    
    dat_cc <- tidyr::drop_na(dat, dplyr::all_of(c("outcome_num", pred_cols)))
    if (nrow(dat_cc) < 4 || length(unique(dat_cc$outcome_num)) < 2) {
      return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
    }
    
    # --- initial model fit
    frm <- stats::as.formula(paste("outcome_num ~", paste(pred_cols, collapse = " + ")))
    fit <- tryCatch(stats::lm(frm, data = dat_cc), error = function(e) NULL)
    if (is.null(fit)) {
      return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
    }
    
    # --- detect and remove collinear variables
    colin_tbl <- tryCatch(detect_collinear_vars(fit), error = function(e) NULL)
    if (!is.null(colin_tbl) && nrow(colin_tbl) > 0) {
      # join PC scores if available
      pc_tbl <- tryCatch(rows |> dplyr::select(variable, pc_neg_logp), error = function(e) NULL)
      
      drops <- purrr::map2_chr(colin_tbl$var1, colin_tbl$var2, function(v1, v2) {
        pc1 <- pc_tbl |> dplyr::filter(variable == v1) |> dplyr::pull(pc_neg_logp)
        pc2 <- pc_tbl |> dplyr::filter(variable == v2) |> dplyr::pull(pc_neg_logp)
        
        # fallback: if both missing, use cardinality
        if (length(pc1) == 0 && length(pc2) == 0) {
          n1 <- if (v1 %in% names(dat_cc)) dplyr::n_distinct(dat_cc[[v1]], na.rm = TRUE) else 1
          n2 <- if (v2 %in% names(dat_cc)) dplyr::n_distinct(dat_cc[[v2]], na.rm = TRUE) else 1
          return(ifelse(n1 >= n2, v1, v2))
        }
        
        # if only one PC value exists, keep that one
        if (is.na(pc1) && !is.na(pc2)) return(v1)
        if (!is.na(pc1) && is.na(pc2)) return(v2)
        
        # both exist — drop the weaker PC association
        ifelse(pc1 >= pc2, v2, v1)
      }) |> unique()
      
      if (length(drops) > 0) {
        message("Dropped collinear (weaker PC) variable(s): ", paste(drops, collapse = ", "))
      }
      
      pred_cols <- setdiff(pred_cols, drops)
      if (length(pred_cols) < 2) {
        return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
      }
      
      frm <- stats::as.formula(paste("outcome_num ~", paste(pred_cols, collapse = " + ")))
      fit <- tryCatch(stats::lm(frm, data = dat_cc), error = function(e) NULL)
      if (is.null(fit)) {
        return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
      }
    }
    
    
    # --- compute VIFs
    v <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (is.null(v)) {
      return(tibble::tibble(variable = covar_cols, vif = NA_real_, gvif_adj = NA_real_))
    }
    
    # --- tidy VIF table
    if (is.null(dim(v))) {
      tibble::tibble(variable = names(v), vif = as.numeric(v), gvif_adj = NA_real_)
    } else {
      v <- as.data.frame(v)
      v$variable <- rownames(v)
      v$gvif_adj <- v$GVIF^(1 / (2 * v$Df))
      tibble::tibble(variable = v$variable, vif = NA_real_, gvif_adj = v$gvif_adj)
    }
  }
  
  
  
  # Tidy per-group processing
  res <- purrr::pmap_dfr(groups_df, function(...) {
    gr <- list(...)
    # subset samples for this group
    idx <- rep(TRUE, nrow(df))
    for (i in seq_along(group_vars)) idx <- idx & (df[[group_vars[i]]] == gr[[i]])
    se_sub <- se[, idx]
    df_sub <- df[idx, , drop = FALSE]
    if (ncol(se_sub) < 3) return(NULL)
    
    pcs <- compute_pcs_for_subset(se_sub)
    if (is.null(pcs)) return(NULL)
    
    # Candidate columns
    cand_cols <- setdiff(colnames(df_sub), exclude_cols)
    if (!is.null(candidate_cols)) cand_cols <- intersect(cand_cols, candidate_cols)
    cand_cols <- cand_cols[!grepl("^PC[0-9]+$", cand_cols)]
    
    # outcome vectors
    outcome_vec <- df_sub[[outcome]]
    outcome_fac <- if (!is.null(outcome_vec) && !is.numeric(outcome_vec)) droplevels(as.factor(outcome_vec)) else NULL
    
    # Evaluate candidates (tidy)
    rows <- purrr::map_dfr(cand_cols, function(cn) {
      v <- df_sub[[cn]]
      if (!is.numeric(v)) {
        v_num <- suppressWarnings(as.numeric(v))
        if (sum(!is.na(v_num)) >= (min_non_na_prop * length(v))) v <- v_num
      }
      
      is_num <- is.numeric(v)
      prop_non_na <- mean(!is.na(v))
      
      if (is_num) {
        sdv <- stats::sd(v, na.rm = TRUE)
        feasible <- is.finite(sdv) && sdv > 0
        n_levels <- NA_integer_
      } else {
        fac <- as.factor(v)
        tbl <- table(fac, useNA = "no")
        n_levels <- length(tbl)
        feasible <- (n_levels >= 2) && (n_levels <= max_levels) && (length(tbl) > 0) && all(tbl >= min_per_level)
        feasible <- isTRUE(feasible)
      }
      
      if (!feasible) {
        return(tibble::tibble(
          !!!gr,
          variable = cn,
          type = if (is_num) "numeric" else "categorical",
          prop_non_na = prop_non_na,
          n_levels = n_levels,
          pc_min_p = NA_real_,
          pc_neg_logp = NA_real_,
          outcome_p = NA_real_,
          feasible = FALSE
        ))
      }
      
      p_pc <- assoc_min_p(v, pcs)
      score <- if (is.finite(p_pc) && !is.na(p_pc)) -log10(p_pc) else NA_real_
      p_out <- safe_p_out(v, outcome_vec, outcome_fac)
      
      tibble::tibble(
        !!!gr,
        variable = cn,
        type = if (is_num) "numeric" else "categorical",
        prop_non_na = prop_non_na,
        n_levels = n_levels,
        pc_min_p = p_pc,
        pc_neg_logp = score,
        outcome_p = p_out,
        feasible = TRUE
      )
    })
    
    if (nrow(rows) == 0) return(NULL)
    
    # VIF: compute on feasible covariates only
    covar_cols_vif <- rows |> dplyr::filter(feasible) |> dplyr::pull(variable) |> unique()
    vif_tbl <- compute_vif_tbl(df_sub, covar_cols_vif, outcome_name = outcome)
    rows |>
      dplyr::left_join(vif_tbl, by = "variable") |>
      dplyr::mutate(
        pc_recommended = pc_neg_logp >= -log10(0.05),
        vif_ok = ifelse(!is.na(vif), vif < 5, ifelse(!is.na(gvif_adj), gvif_adj < 5, NA)),
        include = feasible & pc_recommended & vif_ok & outcome_p > 0.05
      )
  })
  
  if (nrow(res) == 0) return(tibble::tibble())
  
  res |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_vars)),
                   dplyr::desc(pc_recommended),
                   dplyr::desc(pc_neg_logp))
}


# --- Detect Colinear ----------
detect_collinear_vars <- function(fit) {
  ali <- tryCatch(alias(fit)$Complete, error = function(e) NULL)
  if (is.null(ali) || nrow(ali) == 0) {
    return(tibble::tibble(var1 = character(), var2 = character()))
  }
  
  row_names <- rownames(ali)
  col_names <- colnames(ali)
  nz <- which(ali != 0, arr.ind = TRUE)
  if (nrow(nz) == 0) {
    return(tibble::tibble(var1 = character(), var2 = character()))
  }
  
  # helper: remove dummy suffixes to recover original variable name
  get_base <- function(nm) {
    nm <- sub("[:\\.].*", "", nm)          # cut off at first ':' or '.' (factor dummy)
    nm <- sub("[0-9]+$", "", nm)           # remove trailing digits (familyF12 → family)
    trimws(nm)
  }
  
  tibble::tibble(
    aliased_var = row_names[nz[, 1]],
    depends_on  = col_names[nz[, 2]],
    var1 = vapply(row_names[nz[, 1]], get_base, character(1)),
    var2 = vapply(col_names[nz[, 2]], get_base, character(1))
  ) |>
    dplyr::filter(var1 != var2) |>
    dplyr::distinct(var1, var2)
}


# a=suggest_covariates(gene_se,candidate_cols = colnames(gene_se |> pivot_sample())[1:25][!colnames(gene_se |> pivot_sample())[1:25] %in% "dev_stage"],max_levels = 8,log_trans = T,top_var_features = 3000,pc_k = 10,min_per_level = 2)