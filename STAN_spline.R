library("splines")
library("rstan")
require("plot3D")
require(mgcv)
require(visreg)
require(dplyr)
require(forcats)
require(magrittr)


uv <- read.csv("binom_reg_dataset.csv")

N <- nrow(uv)
x1 <- uv$Z
x2 <- scale(uv$STELLAR_MASS)
x3 <-  uv$WHAN_CLASS %>% fct_explicit_na("lineless")

y <- uv$LOGIT_CLASS.1.UVUP.0.UVWEAK.

X <- model.matrix(~ 1 + x1 + I(x1^2) + x2 + I(x2^2) )
K <- ncol(X)
re <- as.numeric(x3)
Nre <- length(unique(re))

# Generated quantities
grid1 <- seq(min(x1),max(x1),length.out = 100)
grid2 <- seq(min(x2),max(x2),length.out = 100)
grid <- expand.grid(grid1,grid2)
xx1 <- grid[,1]
xx2 <- grid[,2]
XX <-  model.matrix(~ 1 + xx1 + I(xx1^2) + xx2 + I(xx2^2) )
M <- nrow(XX)

stan_data  <- list(y = y,
                   X = X,
                   N = N,
                   K = K, 
                   re = re,
                   Nre = Nre,
                   M = M, 
                   XX = XX
                  )



stan_model = "
 data {
  int<lower=1>N;
  int<lower=1>M;
  int<lower=1>Nre;
  int<lower=1> K;
  int re[N];
  matrix[N,K] X; 
  matrix[M,K] XX; 
  int<lower=0, upper=1> y[N];
}

parameters {
    matrix[K,Nre] beta;       // 25 betas!
    real<lower=0> sigma;    // Shared hyperpriors
    real mu;                // Shared hyperpriors

}


model {
   vector[N] pi;
    for (i in 1:N) {
      pi[i] = beta[1,re[i]]*X[i,1] + beta[2,re[i]]*X[i,2] + 
      beta[3,re[i]]*X[i,3] + beta[4,re[i]]*X[i,4] + beta[5,re[i]]*X[i,5];
       }

 // shared hyperpriors
     sigma ~ gamma(0.001, 0.001);
     mu ~ normal(0, 100);

// priors and likelihood
    for (i in 1:K) {
    for (j in 1:Nre) beta[i,j] ~ normal(mu, sigma);
    }

  y ~ bernoulli_logit(pi);
}

generated quantities{
vector[M] pi_1;
vector[M] eta_1;
for(j in 1:M){

eta_1[j] = beta[1,1]*XX[j,1] + beta[2,1]*XX[j,2] + 
      beta[3,1]*XX[j,3] + beta[4,1]*XX[j,4] + beta[5,1]*XX[j,5];
pi_1[j] = inv_logit(eta_1[j]);
}


}

"

fit <- stan(model_code = stan_model,
              data = stan_data,
              seed = 42,
              chains = 3,
              iter = 5000,
              cores = 3,
              warmup = 1500,
               control = list(max_treedepth = 20,
                              adapt_delta=0.99))




pi_1 <- rstan::extract(fit,pars ="pi_1")



p1_quant <- apply(pi_1$pi_1, 2, quantile, c(0.05,0.5, 0.95))

pi_1_mean <- p1_quant[2,]



persp3D(x=grid1, y = sd(uv$STELLAR_MASS)*grid2 + mean(uv$STELLAR_MASS), z = matrix(pi_1_mean,nrow=100,ncol=100, byrow = TRUE),   
        cex = 1, cex.lab=1.5,type="p",pch = 19,alpha=0.5,
        theta = 25, phi = 15, ticktype = "detailed",col="cyan3",bty = "b2",
        xlab="redshift",
        ylab="Log Mass",
        zlab="UV fraction")


