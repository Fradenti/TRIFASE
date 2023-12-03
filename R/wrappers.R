single_run <- function(X, seed, 
                       L, Linv, Dinv, iKER,
                       tau,
                       K, R,
                       upd_row,
                       upd_col,
                       n_iter,
                       nugget,
                       swipes,
                       estimate_tau, 
                       warm=F, HCdX, HCdtX,
                       threshold){
  n <- nrow(X)
  p <- ncol(X)
  
  set.seed(seed)
  Fm <- matrix(0, n, K)
  Gm <- matrix(0, p, R)
  
  if(!warm){
    # completely random start
    for(i in 1:nrow(Fm)) Fm[i, sample(1:ncol(Fm), 1)] <- 1
    for(i in 1:nrow(Gm)) Gm[i, sample(1:ncol(Gm), 1)] <- 1
    # start with hclust
  }else{
    init_cl_row <- stats::cutree(HCdX,k = K)
    init_cl_col <- stats::cutree(HCdtX,k = R)
    Fm <- model.matrix(~ factor(init_cl_row)-1)
    Gm <- model.matrix(~ factor(init_cl_col)-1)
    u <-  sample(1:nrow(Fm),1)
    v <-  sample(1:nrow(Gm),1)
    indF <- sample(1:nrow(Fm),u,replace = FALSE)
    indG <- sample(1:nrow(Gm),v,replace = FALSE)
    for(i in indF) {
      Fm[i,] <- 0
      Fm[i, sample(1:ncol(Fm), 1)] <- 1
    }
    for(i in indG) {
      Gm[i,] <- 0
      Gm[i, sample(1:ncol(Gm), 1)] <- 1
    }
    Fm <- matrix(0, n, K)
    for(i in 1:nrow(Fm)) Fm[i, sample(1:ncol(Fm), 1)] <- 1
  }
  
  
  
  t1 <- Sys.time()
  res <- AlgoAll_cpp2(X = X,
                      Fm = Fm,
                      Gm = Gm,
                      L = L,
                      Linv = Linv,
                      Dinv = Dinv,
                      iKernel = iKER,
                      current_tau = tau,
                      K = K, R = R,
                      upd_row = upd_row,
                      upd_col = upd_col,
                      niter = n_iter,
                      estimate_tau = estimate_tau,
                      nugget = nugget,
                      swipes = swipes,
                      threshold = threshold)
  t2 <- Sys.time()
  f <- apply(res$F, 1, function(x) which(x == 1))
  g <- apply(res$G, 1, function(x) which(x == 1))
  
  return( list(f = f, g = g, Fm = res$F, Gm = res$G, 
               mu = res$mu, loss = res$loss,
               elapsed_time = as.numeric(t2-t1,units="secs"),
               tau2 = res$tau_est,
               ind_minloss = res$ind_minloss,
               niter = length(res$loss)) )
}

#' TRIFASE multi run
#' @importFrom Rfast Dist
#' @export
#' @examples
#' ...
#' 
multi_runs <- function(RUNS=50,CORES = 1,
                       X, 
                      distX = NULL,
                      disttX = NULL,
                      S, 
                      K, 
                      R, 
                      tau, 
                      phi, 
                      kernel = c("exponential","none"),
                      n_iter = 1000,
                      seed = 213,
                      nugget = 1e-8,
                      warm = FALSE,
                      swipes=3,
                      estimate_tau = TRUE,
                      threshold = 0.001,
                      upd_row = c("C","S"),
                      upd_col = c("S","C","A")){
  
  kernel <- match.arg(kernel)
  upd_row <- match.arg(upd_row)
  upd_col <- match.arg(upd_col)
  if(estimate_tau & kernel != "exponential"){
    stop("Nope! You cannot estimate tau without exp kernel")
  }
  
  # initialization
  n <- nrow(X)
  p <- ncol(X)
  if(kernel == "exponential"){
    KER = exp(-as.matrix(Rfast::Dist(S))/phi)
    D <- tau*KER
    iKER = solve(KER)
    cat("Common kernel has been computed!\n")
  }else{
    # placeholder, to be changed:   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    KER <- diag(p)
    D <- KER
    iKER = solve(KER)
  }
  if(warm){
    HCdX = hclust(distX)
    HCdtX = hclust(disttX)
  }
  L    <- t(chol(D))   # lower triangular
  Linv <- solve(L)
  Dinv <- t(Linv) %*% Linv
  cat("Starting parallel runs!\n")
  
  RESULTS <- parallel::mclapply(1:RUNS,function(y){
                single_run(X = X,
                           seed = seed*y, 
                           L = L,
                           Linv = Linv, 
                           Dinv = Dinv, 
                           iKER = iKER,
                           tau = tau,
                           K = K, R = R,
                           upd_row = upd_row,
                           upd_col = upd_col,
                           n_iter = n_iter,
                           nugget = nugget,
                           swipes =swipes, 
                           estimate_tau=estimate_tau,
                           warm = warm, 
                           HCdX = HCdX, 
                           HCdtX = HCdtX,
                           threshold = threshold)},
                mc.cores = CORES)
  return(RESULTS)
}


#' @export
extract_best_min_loss <- function(RESULTS){
  ind <- which.min(unlist(lapply(RESULTS,function(x) min(x$loss))))
  RESULTS[[ind]]
}
