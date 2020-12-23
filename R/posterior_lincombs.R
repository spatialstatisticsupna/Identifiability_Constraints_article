## Posterior distribution of the global intercept ##
lc1 = inla.make.lincomb(Predictor = rep(1/(S*T),S*T))
names(lc1) <- "intercept"


## Posterior distribution of the spatial pattern ##
M1 <- matrix(-1/(S*T),S,S*T)
for (i in 1:S) {
	M1[i,seq(0,S*(T-1),by=S)+i] <- 1/(S*T)*(S-1)
}
lc2 = inla.make.lincombs(Predictor = M1)
names(lc2) <- paste("spatial.",as.character(seq(1:S)),sep="")

                       
## Posterior distribution of the temporal pattern ##
M2 <- matrix(-1/(S*T),T,S*T)
for (i in 1:T) {
	M2[i,seq(1,S)+S*(i-1)] <- 1/(S*T)*(T-1)
}
lc3 = inla.make.lincombs(Predictor = M2)
names(lc3) <- paste("temporal.",as.character(seq(1:T)),sep="")


## Posterior distribution of the spatio-temporal pattern ##
M3 <- matrix(1/(S*T),S*T,S*T)
k=1
for (j in 1:T) {
	for (i in 1:S) {
		M3[k,seq(0,S*(T-1),by=S)+i] <- 1/(S*T)*(1-S)
		M3[k,seq(1,S)+S*(j-1)] <- 1/(S*T)*(1-T)
		M3[k,k] <- (S*T-S-T+1)/(S*T)
		k=k+1
	}
}
lc4 = inla.make.lincombs(Predictor = M3)
names(lc4) <- paste("spatio.temporal.",as.character(seq(1:(S*T))),sep="")

all.lc <- c(lc1,lc2,lc3,lc4)