##########################################################################
#### INLA code to fit spatio-temporal CAR models (Goicoa et al, 2018) ####
##########################################################################
rm(list=ls())
library(sf)

## Download and install the R-INLA package (stable or testing version)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


####################################################
##  Load the original data              				  ##
##	- Spanish female breast cancer mortality data ##
####################################################
load("../data/BreastCancer_ESP.Rdata")
head(Data)

S <- length(unique(Data$Area))
T <- length(unique(Data$Year))

t.from <- min(Data$Year)
t.to <- max(Data$Year)


###################################################################################
## Important: Data must be ordered according to the Kronecker product given for  ##
##		        the structure matrix of the space-time interaction random effect   ##
###################################################################################
Data <- Data[order(Data$Year,Data$Area),]

Data.INLA <- data.frame(O=Data$Counts,
                        E=Data$Expected,
                        ID.area=rep(1:S,T),
                        ID.year=rep(1:T,each=S),
                        ID.area.year=seq(1,T*S))


#####################################################
##  Define the spatial structure matrix of a LCAR  ##
#####################################################

## Load the spatial neighbourhood matrix ##
g <- inla.read.graph("../data/Esp_prov_nb.graph")
Q.xi <- matrix(0, g$n, g$n)
for (i in 1:g$n){
  Q.xi[i,i]=g$nnbs[[i]]
  Q.xi[i,g$nbs[[i]]]=-1
}
Q.xi <- as(Q.xi,"Matrix")
Q.Leroux <- diag(S)-Q.xi


#########################################################
##  Define the temporal structure matrix of a RW1/RW2  ##
#########################################################
D1 <- diff(diag(T),differences=1)
Q.gammaRW1 <- as(t(D1)%*%D1,"Matrix")

D2 <- diff(diag(T),differences=2)
Q.gammaRW2 <- as(t(D2)%*%D2,"Matrix")


######################################################################################
##  Define appropriate hyperprior distributions using the "expression()" function   ##
##	- Unif(0,Inf) for standard deviations							                              ##
##	- Unif(0,1) for the spatial smoothing parameter					                        ##
######################################################################################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"


## Compute posterior patterns distributions ##
compute.patterns <- TRUE  ## Set compute.patterns=FALSE if posterior patterns are not required

if(compute.patterns){
  source("../R/posterior_lincombs.R")
}else{
  all.lc <- NULL
}


##########################################################
## Define the 'formula' objects and fit the INLA models ##
##########################################################
strategy <- "gaussian"  ## Set strategy="simplified.laplace" for more accurate results

## Type I interaction and RW1 prior for time ##
f.TypeI.RW1 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw1", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="iid", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif)))

TypeI.RW1 <- inla(f.TypeI.RW1, family="poisson", data=Data.INLA, E=E,
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  lincomb=all.lc,
                  control.inla=list(strategy=strategy))


## Type I interaction and RW2 prior for time ##
f.TypeI.RW2 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw2", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="iid", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=matrix(rep(1:T,S),1,S*T),e=0))

TypeI.RW2 <- inla(f.TypeI.RW2, family="poisson", data=Data.INLA, E=E,
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  lincomb=all.lc,
                  control.inla=list(strategy=strategy))


## Type II interaction and RW1 prior for time ##
R <- kronecker(Q.gammaRW1,diag(S))
r.def <- S
A.constr <- kronecker(matrix(1,1,T),diag(S))
A.constr <- A.constr[-1,]

f.TypeII.RW1 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw1", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
                 constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=A.constr, e=rep(0,S-1)))

TypeII.RW1 <- inla(f.TypeII.RW1, family="poisson", data=Data.INLA, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   lincomb=all.lc,
                   control.inla=list(strategy=strategy))


## Type II interaction and RW2 prior for time ##
R <- kronecker(Q.gammaRW2,diag(S))
r.def <- 2*S
A.constr <- kronecker(matrix(1,1,T),diag(S))
A.constr <- A.constr[-1,]

f.TypeII.RW2 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw2", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=matrix(1:T,1,T),e=0)) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
                 constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=A.constr, e=rep(0,S-1)))

TypeII.RW2 <- inla(f.TypeII.RW2, family="poisson", data=Data.INLA, E=E,
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  lincomb=all.lc,
                  control.inla=list(strategy=strategy))

## Type III interaction and RW1 prior for time ##
R <- kronecker(diag(T),Q.xi)
r.def <- T
A.constr <- kronecker(diag(T),matrix(1,1,S))
A.constr <- A.constr[-1,]

f.TypeIII.RW1 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw1", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
                 constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=A.constr, e=rep(0,T-1)))

TypeIII.RW1 <- inla(f.TypeIII.RW1, family="poisson", data=Data.INLA, E=E,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    lincomb=all.lc,
                    control.inla=list(strategy=strategy))


## Type III interaction and RW2 prior for time ##
R <- kronecker(diag(T),Q.xi)
r.def <- T
A.constr <- kronecker(diag(T),matrix(1,1,S))
A.constr <- A.constr[-1,]

f.TypeIII.RW2 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw1", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
                 constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=A.constr, e=rep(0,T-1)))

TypeIII.RW2 <- inla(f.TypeIII.RW2, family="poisson", data=Data.INLA, E=E,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    lincomb=all.lc,
                    control.inla=list(strategy=strategy))


## Type IV interaction and RW1 prior for time ##
R <- kronecker(Q.gammaRW1,Q.xi)
r.def <- S+T-1
A1 <- kronecker(matrix(1,1,T),diag(S))
A2 <- kronecker(diag(T),matrix(1,1,S))
A.constr <- rbind(A1[-1,],A2[-1,])

f.TypeIV.RW1 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw1", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
                 constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=A.constr, e=rep(0,S+T-2)))

TypeIV.RW1 <- inla(f.TypeIV.RW1, family="poisson", data=Data.INLA, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   lincomb=all.lc,
                   control.inla=list(strategy=strategy))


## Type IV interaction and RW2 prior for time ##
R <- kronecker(Q.gammaRW2,Q.xi)
r.def <- 2*S+T-2
A1 <- kronecker(matrix(1,1,T),diag(S))
A2 <- kronecker(diag(T),matrix(1,1,S))
A.constr <- rbind(A1,A2)
A.constr <- rbind(A1[-1,],A2[-1,])

f.TypeIV.RW2 <- O ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
               f(ID.year, model="rw2", constr=TRUE,
                 hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
                 constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                 extraconstr=list(A=A.constr, e=rep(0,S+T-2)))

TypeIV.RW2 <- inla(f.TypeIV.RW2, family="poisson", data=Data.INLA, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   lincomb=all.lc,
                   control.inla=list(strategy=strategy))


## Save INLA models ##
MODELS <- list(TypeI.RW1=TypeI.RW1,
               TypeII.RW1=TypeII.RW1,
               TypeIII.RW1=TypeIII.RW1,
               TypeIV.RW1=TypeIV.RW1,
               TypeI.RW2=TypeI.RW2,
               TypeII.RW2=TypeII.RW2,
               TypeIII.RW2=TypeIII.RW2,
               TypeIV.RW2=TypeIV.RW2)

save(MODELS, file="CARmodels_INLA.Rdata")
