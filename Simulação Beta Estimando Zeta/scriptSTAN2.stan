// data block
data{
  int<lower=1> n;
  vector[n] y;
  vector[n] x1;
  vector[n] x2;  
  real m0;
  real m1;
  real m2;
  real<lower=0> v0;
  real<lower=0> v1;
  real<lower=0> v2;
  vector[n] w1;
  vector[n] w2;  
  real mg0;
  real mg1;
  real mg2;
  real<lower=0> vg0;
  real<lower=0> vg1;
  real<lower=0> vg2;
}

// parameters block
parameters{
  real beta0;
  real beta1;
  real beta2;
  real gama0;
  real gama1;
  real gama2;  
}

// Bloco de parametros sem priori
transformed parameters{
  vector[n] theta;
  vector[n] zeta;
  real auxt;
  real auxz;
  for(i in 1:n){
   auxt = beta0 + beta1*x1[i] + beta2*x2[i];
   theta[i] = exp(auxt)/(1+exp(auxt));
  }
  for(i in 1:n){
   auxz = gama0 + gama1*w1[i] + gama2*w2[i];
   zeta[i] = exp(auxz);
  }
}

// model block
model{
  // verossimilhança
  for(i in 1:n){
    y[i] ~ beta(zeta[i]*theta[i],zeta[i]*(1-theta[i]));
  }
  // distribuições a priori  
  beta0 ~ normal(m0,sqrt(v0));
  beta1 ~ normal(m1,sqrt(v1));
  beta2 ~ normal(m2,sqrt(v2));
  gama0 ~ normal(mg0,sqrt(vg0));
  gama1 ~ normal(mg1,sqrt(vg1));
  gama2 ~ normal(mg2,sqrt(vg2));
}
// deixe a linha final abaixo vazia (caso contrário o Stan reclama).
