require(ggplot2)
require(mgcv)
require(xkcd)
# read data

uv <- read.csv("logit_dataset.csv")

set.seed(1056)                           # set seed to replicate example
nobs = nrow(uv)                             # number of obsservations
x1 <- uv$redshift                    # random uniform variable
by <- uv$logit_class.1.uvup.2.uvweak.      # create y as adjusted random bernoulli variate


logitmod <- data.frame(x1,by)


logitgam1 <- gam(by ~ s(x1,bs="cr",k=12),family= binomial(link="logit"),method="REML")

# Bin data for visualization
binx <- 0.05
t.breaks <- cut(x1, seq(min(x1),max(x1), by=binx))
means <- tapply(by, t.breaks, mean)
semean <- function(x) sd(x)/sqrt(length(x))
means.se <- tapply(by, t.breaks, semean)

# Plot
Preds_nzero <- predict(logitgam1,type="link",se=T,unconditional=T) 
fit.link    <- Preds_nzero$fit 
se          <- Preds_nzero$se  
CI.L        <- fit.link-2*se 
CI.R        <- fit.link+2*se 
CI          <- cbind(fit.link,CI.L,CI.R) 
CI          <- exp(CI)/(1+exp(CI)) # The first column correponds to the estimated probability of being non-zero.
colnames(CI) <- c("Predictions","CI_L","CI_R")
CI <- as.data.frame(CI)

gbin <- data.frame(x=seq(binx+min(x1),max(x1), by=binx),y=means)

ggplot(logitmod,aes(x=x1,y=by))+ 
 # geom_point(colour="orange",size=1.25,alpha=0.85,position = position_jitter (h = 0.075))+
  geom_point(aes(x=x,y=y),size=3,data=gbin,colour="cyan3")+
  geom_errorbar(data=gbin,aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),
                colour="cyan3",width=0.01)+
  geom_ribbon(data=CI,aes(ymin=CI_L, ymax=CI_R), alpha=0.45, fill=c("purple"),show.legend=FALSE) +
  geom_line(data=CI,aes(x=x1,y=Predictions),colour="blue3",linetype="dashed",size=1,show.legend=FALSE)+
  theme_xkcd()+xlab("Redshift") + ylab("UV/Red-sequence")+
  theme(axis.title=element_text(size=25)) +
  theme(panel.background = element_rect(color = "black", fill = "gray85") )


require(randomForest)
rby <- as.factor(by)
Rmodel <- randomForest(rby ~ x1,ntree=100)
Rpred <- predict(Rmodel,type = "prob")


ggplot(logitmod,aes(x=x1,y=by))+ 
  geom_point(colour="orange",size=1.25,alpha=0.85,position = position_jitter (h = 0.075))+
  geom_point(aes(x=x,y=y),size=3,data=gbin,colour="cyan3")+
  geom_errorbar(data=gbin,aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),
                colour="cyan3",width=0.01)+
  geom_line(aes(x=x1,y=Rpred[,2]),colour="purple",linetype="dashed",size=1,show.legend=FALSE)+
  theme_xkcd()+coord_cartesian(ylim=c(0,1))+xlab(expression(log~x[mol]))+ylab("Probability of Star Formation")+
  theme(axis.title=element_text(size=25)) +
  theme(panel.background = element_rect(color = "black", fill = "gray85") )






library(neuralnet)
nnd <- data.frame(x1,by)
nn <- neuralnet(by ~ x1,data=nnd,hidden=c(5,3),linear.output = T)

Loessmodel <- loess(by ~ x1,degree = 5)
Lpred <- predict(Loessmodel)


ggplot(logitmod,aes(x=x1,y=by))+ 
  geom_point(colour="orange",size=1.25,alpha=0.85,position = position_jitter (h = 0.075))+
  geom_point(aes(x=x,y=y),size=3,data=gbin,colour="cyan3")+
  geom_errorbar(data=gbin,aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),
                colour="cyan3",width=0.01)+
  geom_line(aes(x=x1,y=Lpred),colour="purple",linetype="dashed",size=1,show.legend=FALSE)+
  theme_xkcd()+coord_cartesian(ylim=c(0,1))+xlab(expression(log~x[mol]))+ylab("Probability of Star Formation")+
  theme(axis.title=element_text(size=25)) +
  theme(panel.background = element_rect(color = "black", fill = "gray85") )
