run_fenr <- function(
    selected_genes,
    universe_genes,
    term_data=NULL,
    db = c("GO", "KEGG", "Reactome", "BioPlanet", "WikiPathways"),
    species = NULL,
    feature_col = "gene_symbol",
    go_url      = "http://current.geneontology.org/annotations"
) {
  options(GO_ANNOTATION_URL = go_url)
  stopifnot(requireNamespace("fenr", quietly = TRUE))
  
  db <- match.arg(db)
  
  fetch_fun <- switch(
    db,
    GO           = fenr::fetch_go,
    KEGG         = fenr::fetch_kegg,
    Reactome     = fenr::fetch_reactome,
    BioPlanet    = fenr::fetch_bioplanet,
    WikiPathways = fenr::fetch_wikipathways
  )
  
  if (db == "GO" && is.null(species)) {
    stop("Please specify a species designation for GO. Use `fetch_go_species()` to see options.")
  }
  
  if(is.null(term_data)){
    term_data <- if (db == "GO") fetch_fun(species = species) else fetch_fun()
  } else {
    term_data
  }
  
  terms_obj <- fenr::prepare_for_enrichment(
    terms        = term_data$terms,
    mapping      = term_data$mapping,
    all_features = universe_genes,
    feature_name = feature_col
  )
  
  fenr::functional_enrichment(
    feat_all  = universe_genes,
    feat_sel  = selected_genes,
    term_data = terms_obj
  )
}

run_fgsea <- function(
    ranked_genes,
    db          = c("GO", "KEGG", "Reactome", "BioPlanet", "WikiPathways"),
    species     = NULL,
    term_data   = NULL,
    feature_col = "gene_symbol",
    go_url      = "https://ftp.ebi.ac.uk/pub/databases/GO/goa",
    min_size    = 10,
    max_size    = 500,
    eps         = 0,
    nperm_seed  = 42
) {
  stopifnot(requireNamespace("fenr",  quietly = TRUE))
  stopifnot(requireNamespace("fgsea", quietly = TRUE))
  
  db <- match.arg(db)
  
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  if (is.null(term_data)) {
    options(GO_ANNOTATION_URL = go_url)
    fetch_fun <- switch(
      db,
      GO           = fenr::fetch_go,
      KEGG         = fenr::fetch_kegg,
      Reactome     = fenr::fetch_reactome,
      BioPlanet    = fenr::fetch_bioplanet,
      WikiPathways = fenr::fetch_wikipathways
    )
    if (db == "GO" && is.null(species))
      stop("Please specify a species for GO. Use `fetch_go_species()` to see options.")
    term_data <- if (db == "GO") fetch_fun(species = species) else fetch_fun()
  }
  
  if (!feature_col %in% names(term_data$mapping))
    stop(
      "feature_col '", feature_col, "' not found. ",
      "Available: ", paste(names(term_data$mapping), collapse = ", ")
    )
  
  gene_sets  <- split(term_data$mapping[[feature_col]], term_data$mapping$term_id) |>
    purrr::map(unique)
  term_names <- term_data$terms |>
    dplyr::distinct(term_id, term_name)
  
  set.seed(nperm_seed)
  fgsea::fgsea(
    pathways = gene_sets,
    stats    = ranked_genes,
    minSize  = min_size,
    maxSize  = max_size,
    eps      = eps
  ) |>
    as.data.frame() |>
    dplyr::mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ",")) |>
    dplyr::rename(term_id = pathway) |>
    dplyr::left_join(term_names, by = "term_id") |>
    dplyr::relocate(term_id, term_name)
}
