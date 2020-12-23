##############################################################
#### Analyse the results of previously fitted INLA models ####
##############################################################
rm(list=ls())
library(INLA)
library(sf)

load("CARmodels_INLA.Rdata")


#####################################################
## Model selection criteria and computational time ##
#####################################################
DIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance, ## posterior mean deviance
             p.eff=x$dic$p.eff,                 ## effective number of parameters
             DIC=x$dic$dic,                     ## Deviance Information Criterion
             WAIC=x$waic$waic,                  ## Watanabe-Akaike information criterion
             LS=-sum(log(x$cpo$cpo)),           ## Logarithmic Score (see inla.cpo function)
             Time=x$cpu.used[4])
}

do.call(rbind,lapply(MODELS, DIC))


################################################################
## Model parameter estimates and standard errors              ##
## (see "inla.hyperpar" function to improve the estimations)  ##
################################################################
Model <- MODELS$TypeIV.RW1

par = c("alpha","sigma2_S","lambda_S","sigma2_T","sigma2_ST")

mean.model <- c(summary(Model)$fixed[1],
                inla.emarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area"),
                Model$summary.hyperpar[2,1],
                inla.emarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.year"),
                inla.emarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area.year"))

sd.model <- c(summary(Model)$fixed[2],
              sqrt(inla.emarginal(function(x) 1/(x^2), Model$marginals.hyperpar$"Precision for ID.area")-mean.model[2]^2),
              Model$summary.hyperpar[2,2],
              sqrt(inla.emarginal(function(x) 1/(x^2), Model$marginals.hyperpar$"Precision for ID.year")-mean.model[4]^2),
              sqrt(inla.emarginal(function(x) 1/(x^2), Model$marginals.hyperpar$"Precision for ID.area.year")-mean.model[5]^2))

table <- data.frame(par, mean.model, sd.model)
print(table)


################################################################################################
##  Percentage of explained variability by the spatial, temporal and spatio-temporal patterns ##
##  (only works if posterior distribution of patterns have been computed)                     ##
################################################################################################
alpha <- Model$summary.lincomb.derived$mean[1]
risks.S <- matrix(Model$summary.lincomb.derived$mean[2:(S+1)],S,T)
risks.T <- t(matrix(Model$summary.lincomb.derived$mean[(S+2):(S+T+1)],T,S))
risks.ST <- matrix(Model$summary.lincomb.derived$'mean'[(S+T+2):(S+T+S*T+1)], S, T)

varS <-var(as.vector(risks.S))
varT <- var(as.vector(risks.T))
varST <- var(as.vector(risks.ST))
c(varS,varT,varST)/(varS+varT+varST)


##############################################
##  Graphical representation of the results ##
##############################################
library(RColorBrewer)
library(tmap)


## Load the cartography of Spanish provinces ##
load("../data/BreastCancer_ESP.Rdata")
plot(Carto_ESP$geometry, axes=TRUE)

S <- length(unique(Data$Area))
T <- length(unique(Data$Year))

t.from <- min(Data$Year)
t.to <- max(Data$Year)


## Posterior mean estimates of the spatial pattern ##
Carto_ESP$spatial.risk <- unlist(lapply(Model$marginals.lincomb.derived[2:(S+1)], function(x) inla.emarginal(exp,x)))
                                           
color.pal <- brewer.pal(6,"RdYlGn")[6:1]
values <- c(0.77,0.83,0.91,1,1.1,1.2,1.3)

Spatial.risk <- tm_shape(Carto_ESP) +
  tm_polygons(col="spatial.risk", palette=color.pal, title="RR", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Spatial pattern", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            outer.margins=c(0.02,0.01,0.02,0.01))

tmap_mode("plot")
print(Spatial.risk)


## Posterior exceedence probabilities of the spatial pattern ##
Carto_ESP$spatial.prob <- unlist(lapply(Model$marginals.lincomb.derived[2:(S+1)], function(x){1-inla.pmarginal(0, x)}))

color.pal <- brewer.pal(6,"Blues")[-1]
values <- c(0,0.1,0.2,0.8,0.9,1)

Spatial.prob <- tm_shape(Carto_ESP) +
  tm_polygons(col="spatial.prob", palette=color.pal, title="Probs", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Spatial pattern", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            outer.margins=c(0.02,0.01,0.02,0.01))

tmap_mode("plot")
print(Spatial.prob)


## Temporal pattern of mortality risks ##
temporal <- unlist(lapply(Model$marginals.lincomb.derived[(S+2):(S+T+1)], function(x) inla.emarginal(exp,x)))
aux <- lapply(Model$marginals.lincomb.derived[(S+2):(S+T+1)], function(x) inla.tmarginal(exp,x))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))

x <- 1:T
plot(range(x),c(0.7,1.3),type="n",xlab="",ylab="", xaxt="n", main="Temporal pattern")
axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(temporal)
abline(h=1,lty=2)


## Maps with posterior median estimates of relative risks ##
Risk <- matrix(Model$summary.fitted.values[1:(S*T),"0.5quant"],S,T)
Risk.df <- as.data.frame(Risk)
colnames(Risk.df) <- paste("Risk",seq(t.from,t.to),sep=".")

Risk.df$Area <- Carto_ESP$Area
Carto <- merge(Carto_ESP, Risk.df, by="Area")

color.pal <- brewer.pal(6,"RdYlGn")[6:1]
values <- c(-Inf,0.77,0.91,1,1.10,1.3,Inf)

Map.Risk <- tm_shape(Carto) +
  tm_polygons(col=paste("Risk",seq(t.from,t.to,by=2),sep="."), palette=color.pal, title="Risks", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Breast cancer mortality data", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            panel.labels=seq(t.from,t.to,by=2),
            outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=3, ncol=4)

tmap_mode("plot")
print(Map.Risk)


## Maps with posterior exceedence probabilities of relative risks ##
Prob <- matrix(1-Model$summary.fitted.values[1:(S*T),"1 cdf"],S,T)
Prob.df <- as.data.frame(Prob)
colnames(Prob.df) <- paste("Prob",seq(t.from,t.to),sep=".")

Prob.df$Area <- Carto_ESP$Area
Carto <- merge(Carto_ESP, Prob.df, by="Area")

color.pal <- brewer.pal(6,"Blues")[-1]
values <- c(0,0.1,0.2,0.8,0.9,1)

Map.Prob <- tm_shape(Carto) +
  tm_polygons(col=paste("Prob",seq(t.from,t.to,by=2),sep="."), palette=color.pal, title="Probs", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Breast cancer mortality data", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            panel.labels=seq(t.from,t.to,by=2),
            outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=3, ncol=4)

tmap_mode("plot")
print(Map.Prob)


## Area-specific relative risks ##
q1 <- matrix(Model$summary.fitted.values[1:(S*T),"0.025quant"],S,T)
q2 <- matrix(Model$summary.fitted.values[1:(S*T),"0.975quant"],S,T)
x <- 1:T

plot.areas <- c(8,20,28,31)
plot.SMR <- TRUE

par(mfrow=c(2,2), pty="s")
for(i in plot.areas){
  plot(range(x),c(0.6,1.5),type="n", xlab="", ylab="", xaxt="n", main=Carto_ESP$Name[i])
  axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0)
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(q1[i,], tail(q2[i,], 1), rev(q2[i,]), q1[i,1])
  polygon(X.Vec, Y.Vec, col = "grey", border = NA)
  lines(Risk[i,])
  if(plot.SMR){
    lines(matrix(Model$.args$data$O/Model$.args$data$E,S,T)[i,], col="red", lwd=2)
  }
  abline(h=1,lty=2)
}

