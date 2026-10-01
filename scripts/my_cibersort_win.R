library(future)
library(furrr)
library(e1071)
library(preprocessCore)

CoreAlg <- function(X, y, cores = 3){
  
  
  ########################
  ## X is the data set
  ## y is labels for each row in X
  ########################
  
  
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    
    #if(i==1){nus <- 0.997}
    #if(i==2){nus <- 0.998}
    #if(i==3){nus <- 0.999}
    
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=FALSE)
    model
  }
  
  #Execute In a parallel way the SVM
  if(cores>1){
    if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1)
    else out <-  parallel::mclapply(1:svn_itor, res, mc.cores=cores)
  }
  else out <-  lapply(1:svn_itor, res)
  
  #Initiate two variables with 0
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  ##############################
  ## Here CIBERSORT starts    #
  ##############################
  
  t <- 1
  while(t <= svn_itor) {
    
    #Get the weights with a matrix multiplications between two vectors. I should get just one number (?)
    #This is done multiplying the coefficients (?) and ???
    
    #The support vectors
    #are the points of my dataset that lie closely to the plane that separates categories
    #The problem now is that I don't have any category (discrete variable, e.g., "sport", "cinema") but I ave continuous variable
    mySupportVectors <- out[[t]]$SV
    
    #My coefficients
    myCoefficients <- out[[t]]$coefs
    
    weights = t(myCoefficients) %*% mySupportVectors
    
    #set up weight/relevance on each
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    
    #This multiplies the reference profile for the correspondent weigth
    u <- sweep(X,MARGIN=2,w,'*')
    
    #This does the row sums
    k <- apply(u, 1, sum)
    
    #Don't know
    nusvm[t] <- sqrt((mean((k - y)^2))) #pitagora theorem
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  #print(mn)
  model <- out[[mn]]
  
  #get and normalize coefficients
  
  #############################################
  ## THIS IS THE SECRET OF CIBERSORT
  #############################################
  
  q <- t(model$coefs) %*% model$SV
  
  #############################################
  #############################################
  
  q[which(q<0)]<-0
  
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#' @importFrom stats sd
#'
#' @keywords internal
#'
doPerm <- function(perm, X, Y, cores = 3){
  
  
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, cores = cores)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

# MADE BY STEFANO TO ALLOW PARALLELISM
call_core = function(itor, Y, X, P, pval, CoreAlg, cores = 1){
  ##################################
  ## Analyze the first mixed sample
  ##################################
  
  y <- Y[,itor]
  
  #standardize mixture
  y <- (y - mean(y)) / sd(y)
  
  #run SVR core algorithm
  result <- CoreAlg(X, y, cores = cores)
  
  #get results
  w <- result$w
  mix_r <- result$mix_r
  mix_rmse <- result$mix_rmse
  
  #calculate p-value
  if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
  
  #print output
  c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
  
}


my_cibersort_win <- function(Y, X, perm = 0, QN = TRUE, 
                             workers = parallel::detectCores() - 1,
                             exp_transform = FALSE) {
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  common_genes <- intersect(rownames(X), rownames(Y))
  X <- X[common_genes, , drop = FALSE]
  Y <- Y[common_genes, , drop = FALSE]
  
  P <- perm
  
  if (is.null(exp_transform)) exp_transform <- max(Y) < 50
  if (exp_transform) Y <- 2^Y
  
  if (QN) {
    tmpc <- colnames(Y); tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc; rownames(Y) <- tmpr
  }
  
  Xgns <- rownames(X); Ygns <- rownames(Y)
  Y <- Y[Ygns %in% Xgns, , drop = FALSE]
  X <- X[Xgns %in% rownames(Y), , drop = FALSE]
  
  zero_samples <- colSums(Y) == 0
  if (any(zero_samples))
    warning("Samples with 0 counts removed: ", 
            paste(colnames(Y)[zero_samples], collapse = ", "))
  Y <- Y[, !zero_samples, drop = FALSE]
  
  flat_sd <- matrixStats::colSds(Y) == 0
  if (any(flat_sd))
    warning("Samples with sd = 0 removed: ", 
            paste(colnames(Y)[flat_sd], collapse = ", "))
  Y <- Y[, !flat_sd, drop = FALSE]
  
  X <- (X - mean(X)) / sd(as.vector(X))
  Y_norm <- apply(Y, 2, \(mc) (mc - mean(mc)) / sd(mc))
  
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y)$dist)
  } else {
    nulldist <- NULL
  }
  
  mix <- ncol(Y)
  pval <- 9999
  
  # Works on Windows, Mac, Linux
  plan(multisession, workers = workers)
  
  output <- furrr::future_map(
    1:mix,
    \(itor) call_core(itor, Y, X, P, pval, CoreAlg, cores = 1),
    .options = furrr_options(
      seed = TRUE,
      globals = c("Y", "X", "P", "pval", "CoreAlg", "call_core")
    )
  ) |>
    (\(x) matrix(
      unlist(lapply(x, \(v) v[-1])),  # drop sample name from each vector
      nrow = length(x),
      byrow = TRUE
    ))()
  
  rownames(output) <- colnames(Y)
  colnames(output) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  
  list(proportions = output, mix = Y_norm, signatures = X)
}

update_coldata <- function(se, df) {
  existing_cols <- colnames(tidybulk::pivot_sample(se))
  
  new_coldata <- tidybulk::pivot_sample(se) |>
    dplyr::left_join(
      df |> dplyr::select(-dplyr::any_of(setdiff(existing_cols, ".sample"))),
      by = ".sample"
    ) |>
    tibble::column_to_rownames(".sample")
  
  stopifnot(identical(
    rownames(new_coldata),
    rownames(as.data.frame(SummarizedExperiment::colData(se)))
  ))
  
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(new_coldata)
  se
}
