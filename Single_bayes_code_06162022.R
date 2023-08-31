# This function fit a simple Bayesian linear model for the outcome
bayes_single_model <- function(X, p, y, A, single_mod){
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  D <-  ncol(X)
  
  studydata <- list(
    N = nrow(X), D=D, y = y, 
    X=X,A=A, A_avg = mean(A))
  
  fit_single <- single_mod$sample(
    data = studydata,
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    show_messages = FALSE) 
  
  # Check model divergence and get posterior distributions of parameters 
  diagnostics_df <- as_draws_df(fit_single$sampler_diagnostics())
  div_single <- sum(diagnostics_df[, 'divergent__'])
  
  # Get posterior draws of all parameters
  draws_dt <- data.table(as_draws_df(fit_single$draws()))
  
  #beta_0<- as.matrix(draws_dt[,c("beta_0")]) # estimation of treatment main effect
  beta_1 <- as.matrix(draws_dt%>%select(paste0("beta","[",1:p,"]"))) #estimation of coefficients for treatment by covariate interaction
  m <- as.matrix(draws_dt%>%select(paste0("m","[",1:p,"]"))) #estimation of covariates' main effects 
  #tau <- as.matrix(draws_dt[,c("tau")]) # intercept
  return(data.table(beta_1,
                    m, div_single)) 
}


pred_sng_bayes <- function(xnew=NULL, xmnew=NULL, newA=NULL,A.mean=NULL, coefs = sgbayes.coefs,thresh = 0){
  ##Given patients characteristics, derive optimal treatment for each patients using simple Bayesian linear model
  if(is.null(xnew))  xnew  <- bsim.obj$X
  if(is.null(xmnew)) xmnew <- bsim.obj$Xm
  if(is.null(newA))  newA  <- bsim.obj$A
  
  n <- nrow(xnew)
  p <- ncol(xnew)
  
 # beta_0 <- as.matrix(coefs[,c("beta_0")])
  beta_1 <- as.matrix(coefs%>%select(paste0("beta","[",1:p,"]"))) 
  m <- as.matrix(coefs%>%select(paste0("m","[",1:p,"]"))) #estimation of covariates' main effects 
  #tau <- as.matrix(coefs[,c("tau")]) # int
  
  nsample <-  nrow(beta_1)
  
  eta.distr = contrast.distr = tbi.tmp  <- matrix(rep(0, n*nsample), n, nsample)
  
  newA.mean <- newA-A.mean
  
  for(i in 1:n){
    eta.distr[i,] <- m%*%xnew[i,] + beta_1%*%xnew[i,]*newA.mean[i]#eta.distr (n.test*n.samples): a canonical parameter matrix with a row for each observation(test data point) and a column for each posterior sample
    contrast.distr[i,] <- beta_1%*%xnew[i,]# samples of distribution of tbi (canonical parameter of A=2 - canonical parameter of A=1)
    tbi.tmp[i,] <- (contrast.distr[i,] < thresh)#each row: ? tbi < 0 for the ith patient
  }
  
  tbi <- apply(tbi.tmp, 1, mean)#each row: Pr(tbi < 0) for the ith patient
  
  return(list(eta.distr=eta.distr, contrast.distr =contrast.distr, tbi=tbi))
}

