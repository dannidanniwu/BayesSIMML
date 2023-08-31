# This function fit a simple Bayesian linear model for the outcome
bayes_single_model <- function(X = X, p = p, y = y, A= A, single_mod = single_mod){
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
  
  beta_0<- as.matrix(draws_dt[,c("beta_0")]) # estimation of treatment main effect
  beta_1 <- as.matrix(draws_dt%>%select(paste0("beta","[",1:p,"]"))) #estimation of coefficients for treatment by covariate interaction
  m <- as.matrix(draws_dt%>%select(paste0("m","[",1:p,"]"))) #estimation of covariates' main effects 
  tau <- as.matrix(draws_dt[,c("tau")]) # intercept
  return(data.table(beta_0,beta_1,
                    m, tau, div_single)) 
}


optim_sng_bayes <- function(xnew=NULL, xmnew=NULL, newA=NULL,A.mean=NULL, coefs = sgbayes.coefs,thresh = 0){
  ##Given patients characteristics, derive optimal treatment for each patients using simple Bayesian linear model
  if(is.null(xnew))  xnew  <- bsim.obj$X
  if(is.null(xmnew)) xmnew <- bsim.obj$Xm
  if(is.null(newA))  newA  <- bsim.obj$A
  
  n <- nrow(xnew)
  p <- ncol(xnew)
  
  beta_0 <- as.matrix(coefs[,c("beta_0")])
  beta_1 <- as.matrix(coefs%>%select(paste0("beta","[",1:p,"]"))) 
  m <- as.matrix(coefs%>%select(paste0("m","[",1:p,"]"))) #estimation of covariates' main effects 
  tau <- as.matrix(coefs[,c("tau")]) # int
  
  nsample <-  nrow(beta_0)
  
  eta.distr  <- matrix(rep(0, n*nsample), n, nsample)
  tbi.distr = prob.tbi.tmp <- matrix(rep(0, n*nsample), nsample,n)
  
  tbi.distr <- apply(xnew,1, function(x) beta_0 + beta_1%*%x) # samples of distribution of tbi (canonical parameter of A=2 - canonical parameter of A=1)
  #tbi_distr: a n_sample*n.test matrix: column[i]: posterior draws of beta_0[1]+x[i]*beta_1 for the ith subject
  prob.tbi.tmp <- apply(tbi.distr,2, function(x) x < thresh)
  prob.tbi <- apply(prob.tbi.tmp, 2, mean) #each column: Pr(tbi < 0) for the ith patient
  optTr  <- sapply(prob.tbi, function(x) ifelse(x > 0.5, 2, 1)) #if Pr(tbi < 0) greater than 0.5, recommend CCP
  
  newtrt <- newA-A.mean
  #eta.distr (n.test*n.samples): a canonical paramter matrix with a row for each observation(test data point) and a column for each posterior sample
  for(i in 1:n){
  eta.distr[i,] <- tau + m%*%xnew[i,] + beta_0*newtrt[i] + beta_1%*%xnew[i,]
  }
  
  return(list(tbi.distr = tbi.distr, prob.tbi = prob.tbi, eta.distr = eta.distr, optTr = optTr))
}

