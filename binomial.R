# From: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida
#
# you are kindly asked to include the complete citation if you used this
# material in a publication

# Code 5.25 - Synthetic data from a binomial model in R

set.seed(33559)
source("https://raw.githubusercontent.com/johnbaums/jagstools/master/R/jagsresults.R")

dat <- read.csv("https://raw.githubusercontent.com/mdastro/UV_ETGs/master/Coding/Model/regression_dataset.csv")
nobs <- nrow(dat)

bindata=data.frame(y=dat$number_galaxies_uvup,m=dat$number_galaxies_redsequence,
                   x1 = dat$average_redshift)

# Code 5.26 - Binomial model in R using JAGS

library(R2jags)

X <- model.matrix(~ x1 + I(x1^2), data = bindata)
K <- ncol(X)

model.data <- list(Y = bindata$y,
                   N = nrow(bindata),
                   X = X,
                   K = K,
                   m = bindata$m)

sink("GLOGIT.txt")

cat("
    model{
    # Priors
    # Diffuse normal priors betas
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}
    # Likelihood
    for (i in 1:N){
    Y[i] ~ dbin(p[i],m[i])
    logit(p[i]) <- eta[i]
    eta[i] <- inprod(beta[], X[i,])
    }
# Prediction for new data
 for (j in 1:N){
    etax[j]<-inprod(beta[], X[j,])
    logit(px[j]) <- etax[j]
    Yx[j]~dbern(px[j])
 }
    }
    ",fill = TRUE)

sink()

inits <- function () {list(beta = rnorm(K, 0, 0.1)) }

params <- c("beta","px")

BINL <- jags(data = model.data,
             inits = inits,
             parameters = params,
             model.file = "GLOGIT.txt",
             n.thin = 1,
             n.chains = 3,
             n.burnin = 3000,
             n.iter = 5000)

print(BINL, intervals=c(0.025, 0.975), digits=3)

# Plot
y <- jagsresults(x=BINL, params=c('px'))
x <- bindata$x1
gdata <- data.frame(x = x, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],upr1=y[,"75%"],upr2=y[,"97.5%"])

bindata$frac <- bindata$y/bindata$m
ggplot(bindata,aes(x=x1,y=frac))+
  geom_point(size=2.75,colour="blue3")+
  geom_ribbon(data=gdata,aes(x=x,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.45, fill=c("orange2"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=x,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.35, fill = c("orange"),show.legend=FALSE) +
  geom_line(data=gdata,aes(x=x,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw() +xlab("Redshift") + ylab("UV/Red sequence")



