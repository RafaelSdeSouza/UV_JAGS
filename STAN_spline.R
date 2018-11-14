library("splines")
library("rstan")
require("plot3D")
require(mgcv)
require(visreg)


uv <- read.csv("binom_reg_dataset.csv")

nobs <- nrow(uv)
x1 <- uv$Z
x2 <- uv$STELLAR_MASS
x3 <-  uv$WHAN_CLASS
y <- uv$LOGIT_CLASS.1.UVUP.0.UVWEAK.

X <- model.matrix(~ 1 + x1 + I(x1^2) + x2 + I(x2^2) )
K <- ncol(X)
re <- as.numeric(x3)
Nre <- length(unique(re))

stan_data  <- list(y = y,
                   X = X,
                   x2 = x2,
                   x3 = x3
                  )



stan_model= "

functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
  vector[size(t)] b_spline;
  vector[size(t)] w1 = rep_vector(0, size(t));
  vector[size(t)] w2 = rep_vector(0, size(t));
  if (order==1)
  for (i in 1:size(t))
  b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
  else {
  if (ext_knots[ind] != ext_knots[ind+order-1])
  w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
  (ext_knots[ind+order-1] - ext_knots[ind]);
  if (ext_knots[ind+1] != ext_knots[ind+order])
  w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
  (ext_knots[ind+order] - ext_knots[ind+1]);
  b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
  w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
}
return b_spline;
}
}

data {
  int N;
  int num_knots;
  vector[num_knots] knots;
  int spline_degree; 
  real Y[N];
  real X[N];
}


transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, N] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
  B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[num_knots + spline_degree - 1, N] = 1; 
}



parameters {
row_vector[num_basis] a_raw;
real a0;
real<lower=0> sigma;
real<lower=0> tau;

}



transformed parameters {
  row_vector[num_basis] a;
vector[N] Y_hat;
a[1] = a_raw[1];
for (i in 2:num_basis)
a[i] = a[i-1] + a_raw[i]*tau; 
Y_hat = a0*to_vector(X) + to_vector(a*B);
}


model {
a_raw ~ normal(0, 1);
tau ~ cauchy(0, 1);
sigma ~ cauchy(0, 1);
Y ~ normal(Y_hat, sigma);
}

generated quantities{
vector[N] mu_pred;
// Posterior parameter distribution of the mean

mu_pred = a0*to_vector(X) + to_vector(a*B);

}


"

fit <- stan(model_code = stan_model,
              data = stan_data,
              seed = 42,
              chains = 3,
              iter =1500,
              cores= 3,
              warmup=750)




Y_mean <- extract(fit, "mu_pred")
Y_mean_cred <- apply(Y_mean$mu_pred, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$mu_pred, 2, mean)


fitdat <- data.frame(x = X, y = Y_mean_mean, lwr1 = Y_mean_cred[1,],upr1 = Y_mean_cred[2,])

ggplot(data=gdat,aes(x=x,y=y)) +
  geom_point() +
  geom_line(data=fitdat,aes(x=x,y=y)) +
  geom_ribbon(data=fitdat,aes(x=x,ymin=lwr1, ymax=upr1,y=NULL),fill=c("gray50"),alpha=0.5) +
  theme_bw()

