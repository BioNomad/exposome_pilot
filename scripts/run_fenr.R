run_fenr <- function(
    selected_genes,
    universe_genes,
    db = c("GO", "KEGG", "Reactome", "BioPlanet", "WikiPathways"),
    species = NULL,
    feature_col = "gene_symbol"
) {
  stopifnot(requireNamespace("fenr", quietly = TRUE))
  
  db <- match.arg(db)
  
  # Fetch functional terms + mapping
  fetch_fun <- switch(db,
                      GO = fenr::fetch_go,
                      KEGG = fenr::fetch_kegg,
                      Reactome = fenr::fetch_reactome,
                      BioPlanet = fenr::fetch_bioplanet,
                      WikiPathways = fenr::fetch_wikipathways
  )
  
  # Handle species if required
  if (db == "GO" && is.null(species)) {
    stop("Please specify a species designation for GO. Use `fetch_go_species()` to see options.")
  }
  
  if (db == "GO") {
    term_data <- fetch_fun(species = species)
  } else {
    term_data <- fetch_fun()
  }
  
  # Prepare for enrichment
  terms_obj <- fenr::prepare_for_enrichment(
    terms = term_data$terms,
    mapping = term_data$mapping,
    all_features = universe_genes,
    feature_name = feature_col
  )
  
  # Run enrichment
  enrichment_results <- fenr::functional_enrichment(
    feat_all = universe_genes,
    feat_sel = selected_genes,
    term_data = terms_obj
  )
  
  return(enrichment_results)
}
