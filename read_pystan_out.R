require(magrittr);require(dplyr)

redshift <- read.csv("/Users/rafaeldesouza/Downloads/model_prob.csv")[,2]

malu <- read.csv("/Users/rafaeldesouza/Downloads/fit_results.csv") %>%
  filter(.,grepl("^pnew",parameter))

y <- malu
x <- redshift
gdata <- data.frame(x =redshift, mean = y[,"mean"],lwr1=y[,"X25."],lwr2=y[,"X2.5."],upr1=y[,"X75."],upr2=y[,"X97.5."])





ggplot(gdata,aes(x=redshift,y=mean))+
  geom_ribbon(data=gdata,aes(x=redshift,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.95, fill = c("#fac901"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=redshift,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("#225095"),show.legend=FALSE) +
  geom_line(data=gdata,aes(x=redshift,y=mean),colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw() + xlab("Redshift") + ylab("UV/Red-sequence")