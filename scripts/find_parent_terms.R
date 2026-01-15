find_parent_terms <- function(
    fenr_res,
    threshold=0.7
){
  library(rrvgo)
  
  combined_res <- data.frame()
  
  for(ont in c("BP","CC","MF")){
    message("Calculating Similarity for:",ont)
    simMatrix <- calculateSimMatrix(fenr_res$term_id,
                                    orgdb="org.Hs.eg.db",
                                    ont=ont,
                                    method="Rel")
    
    scores <- setNames(-log10(fenr_res$p_adjust), fenr_res$term_id)
    
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=threshold,
                                    orgdb="org.Hs.eg.db")
    combined_res <- combined_res |> 
      bind_rows(reducedTerms)
  }
  
  return(combined_res)
}
