# Astronomical Data Analysis Summer School
# Bayesian tutorial by Rafael S. de Souza - ELTE, Hungary & COIN
#
# Partial example from Bayesian Models for Astrophysical Data
# by Hilbe, de Souza & Ishida, 2016, Cambridge Univ. Press
#
# Example of Bayesian Bernoulli regression in R using JAGS
# synthetic data
# 1 response (y) and 1 explanatory variable (x1)

require(R2jags)
require(ggplot2)
source("https://raw.githubusercontent.com/johnbaums/jagstools/master/R/jagsresults.R")
uv <-read.csv("logit_dataset.csv")

set.seed(1056)                           # set seed to replicate example
nobs= nrow(uv)                             # number of obsservations
x1 <- uv$redshift                    # random uniform variable
by <- uv$logit_class.1.uvup.2.uvweak.      # create y as adjusted random bernoulli variate


# Prepare data for prediction
M=500
xx = seq(from =  min(x1),
         to =  max(x1),
         length.out = M)


# Construct data dictionary
logitmod <-data.frame(by, x1)
X <- model.matrix(~ 1+x1,
                  data = logitmod)
K <- 3
logit_data <- list(Y  = logitmod$by, # Response variable
                   X  = X,           # Predictors
                   K  = K,           # Number of Predictors including the intercept
                   N  = nobs,        # Sample size
                   xx = xx,
                   M = M
)

# JAGS code
LOGIT<-"model{

# Diffuse normal priors for predictors
for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }

# Likelihood function
for (i in 1:N){
Y[i] ~ dbern(p[i])
logit(p[i]) <- eta[i]
eta[i]  <- beta[1]+beta[2]*X[i,2]+beta[3]*X[i,2]^2
}

# Prediction for new data
for (j in 1:M){
etax[j]<-beta[1]+beta[2]*xx[j]+beta[3]*xx[j]^2
logit(px[j]) <- etax[j]
Yx[j]~dbern(px[j])
}
}"

#A function to generate initial values for mcmc
inits  <- function () {
  list(beta  = rnorm(3, 0, 0.1)  )
}

# define parameters
params <- c("beta","px")

# Fit
jagsfit<- jags(data       = logit_data,
               inits      = inits,
               parameters = params,
               model      = textConnection(LOGIT),
               n.thin     = 1,
               n.chains   = 3,
               n.burnin   = 5000,
               n.iter     = 10000)
#traplot(jagsfit,c("beta"))

# check results
print(jagsfit,intervals=c(0.025, 0.975),justify = "left", digits=2)


# Plot
y <- jagsresults(x=jagsfit, params=c('px'))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],upr1=y[,"75%"],upr2=y[,"97.5%"])

# Bin data for visualization
binx<-0.1
t.breaks <-cut(x1, seq(min(x1),max(x1), by=binx))
means <-tapply(by, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(by, t.breaks, semean)


gbin<-data.frame(x=seq(binx,max(x1), by=binx),y=means)

ggplot(logitmod,aes(x=x1,y=by))+
  geom_point(colour="#dd0100",size=1,alpha=0.45,position = position_jitter (h = 0.025))+
 geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.95, fill = c("#fac901"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("#225095"),show.legend=FALSE) +
  geom_point(aes(x=x,y=y),size=2.75,data=gbin,colour="#dd0100")+
  geom_errorbar(data=gbin,aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.85,
                colour="#dd0100",width=0.005)+
  geom_line(data=gdata,aes(x=xx,y=mean),colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw() + xlab("Redshift") + ylab("UV/Red-sequence")

