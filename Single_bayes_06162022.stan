//no intercept and no main effect of A
data {
  int<lower=0> N;
  int<lower=1> D;
  real A_avg;
  int<lower= 1,upper= 2> A[N];  //treatment assignment
  int<lower=0,upper=1>  y[N];               
  row_vector[D] X[N];       //covariates
}

transformed data {
  vector[N] trt;
  for (i in 1:N)
  trt[i] = A[i] - A_avg;
}

parameters {
 // real tau;            //itercept for logistic model
  unit_vector[D] m; //2->K
 // real beta_0;
  unit_vector[D] beta;      //y specific scale
}
transformed parameters {
  vector[N] yhat;
  
  //for (k in 1:K)
  for (i in 1:N){  
    yhat[i] = X[i]*m + X[i] * beta*trt[i];
  }
  
}
model {
  //priors
  
  beta ~ normal (0,1);
  m ~ normal(0,1);
  //beta_0 ~ normal (0,10); //normal  diffuse prior for beta_0
  //tau ~ student_t(3, 0, 8);
  
      
for (i in 1:N){
    y[i] ~ bernoulli_logit(yhat[i]); 
}
}

