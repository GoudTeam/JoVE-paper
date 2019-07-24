##########################################################################################
##Spatial analysis of lysosomal exocytosis
##By: Hugo LACHUER (hugo.lachuer@curie.fr)
## (R version 3.4.1 (2017-06-30) -- "Single Candle")
## Last update: 24/07/2019
##########################################################################################

##########################################################################################
## 1/ User part:

DatasetPath <- "C:/Users/.../DemoDataset"
PlotPath <- "C:/Users/Bruno/.../Plots"
Ripley <- FALSE	#TRUE --> generate Ripley's K function for each cell, FALSE --> skip this long step
RandomSampling <- TRUE	#TRUE --> generate random re-sampling for original data to give the same statistical weight to each cell, FALSE --> do not make random re-sampling

##########################################################################################
## 2/ Library:

usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("spatstat")
usePackage("DescTools")
usePackage("evmix")

library(spatstat)
library(DescTools)
library(evmix)

##########################################################################################
## 3/ Settings:

setwd(DatasetPath)
parameter <- read.delim("Spherical parameter.txt", header=TRUE)
pattern <- as.data.frame(read.table("Pattern parameter.txt", sep="\t", header=TRUE))
n <- nrow(parameter)
Removed <- 0
META <- as.data.frame(matrix(, ncol=5, nrow=1))

graduation1 <- c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
graduation2 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
vec.expr1 <- parse(text = graduation1)
graduation3 <- c("0", "0.5", "1", "1.5", "2")
graduation4 <- c(0, 0.5, 1, 1.5, 2)
vec.expr2 <- parse(text = graduation3)
SEM_pattern <- sd(pattern[,4])/sqrt(nrow(pattern))

N_sim <- 1000	#Number of Monte-carlo simulation for confidence band on the modulus distribution
Echan <- 1000
NBoot <- 300	#Number of bootstrap
evalpoints <- 1000
N <- 10000 #Number of Monte-Carlo simulation for CSR test based on NND
alpha <- 0.95	#Percentage for error and confidence band

# discrete approximation to the unit disk
w <- owin(c(-1,1), c(-1,1), mask=matrix(TRUE, 1000,1000))
X <- raster.x(w)
Y <- raster.y(w)
cellwindow <- owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))

#Functions
errorbar <- function(x, y, xmax, xmin, h) {
	segments(xmin, y, x1 = xmax, y1 = y, lwd=2)
	segments(xmin, y-h/2, x1 = xmin, y1 = y+h/2, lwd=2)
	segments(xmax, y-h/2, x1 = xmax, y1 = y+h/2, lwd=2)
}

#List of K function
KfunctionList <- list()

##########################################################################################
## 4/ Singe cell spatial analysis:

for(i in 1:n){
	
	##Open file
	setwd(DatasetPath)
	ID <- parameter[i,1]
	Name <- paste0("Results corriged (", ID, ").txt")
	Data <- as.data.frame(read.table(Name, sep="\t", header=TRUE))

	##Polar coordinate
	radius <- parameter[i,5]
	x <- parameter[i,2]
	y <- parameter[i,3]
	Data[,2] <- Data[,2]-x
	Data[,3] <- (Data[,3]-y)*(-1)
	Data2 <- as.data.frame(matrix(0, ncol=4, nrow=nrow(Data)))
	colnames(Data2) <- c("time", "Theta", "A", "A normalized")
	Data2[,1] <- Data[,1]	
	angle <- vector(, nrow(Data))
	for(j in 1:nrow(Data)) {
		Data2[j,3] <- sqrt(Data[j,2]^2 + Data[j,3]^2)
		Data2[j,4] <- Data2[j,3]/radius
		if(Data[j,3] >= 0) {Data2[j,2] <- acos(Data[j,2]/Data2[j,3]) * 180/pi}
		if(Data[j,3] < 0) {Data2[j,2] <- acos(-Data[j,2]/Data2[j,3]) * 180/pi +180}
		if(Data[j,3] >= 0) {angle[j] <- acos(Data[j,2]/Data2[j,3])}
		if(Data[j,3] < 0) {angle[j] <- acos(-Data[j,2]/Data2[j,3]) + pi}
	}
	##Remove cell with A normalized > 1
	for(j in nrow(Data2):1) {
		if(Data2[j,4]>1) {
			Data2 <- Data2[-j,]
			Removed <- Removed + 1
		}
	}
	print(paste0(Removed, "cells removed from ", ID))
	Removed <- 0
	
	setwd(PlotPath)
	
	##Scatter plot (single cell)
	pattern_radius <- as.numeric(as.character(pattern[i,4]))	
	pdf(paste0("Scatter plot (", ID, ").pdf"))
	Coordinates_Normalized <- Data[,(2:3)]/radius
	plot(Coordinates_Normalized, xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), pch=16, col="red", lwd=3, asp=1,  main= paste0("Exocytosis scatter plot (", ID, ") (n=", nrow(Data2), ")"))
	DrawCircle(x = 0, y = 0, r.out = 1, r.in = 0, theta.1 = 0, theta.2 = 2*pi, border = par("fg"), lty = par("lty"), lwd = 3, nv = 1000, plot = TRUE)
	DrawCircle(x = 0, y = 0, r.out = 1, r.in = (1-pattern_radius), theta.1 = 0, theta.2 = 2*pi, border = NA, lty = par("lty"), lwd = 3, nv = 1000,col= rgb(0.06,0.306,0.545,alpha=0.40), plot = TRUE)
	dev.off()
	
	##2D density plot
	lyso <- ppp(Data[,2]/radius,Data[,3]/radius,window=cellwindow)
	lyso.dens <- density(lyso)
	pdf(paste0("Density plot (", ID, ").pdf"))
	plot(lyso.dens, main=paste0("Density plot (", ID, ")"))
	contour(lyso.dens, add=T)
	points(lyso, pch=20)
	dev.off()
	
	##plot Modulus distribution (single cell)
	a <- as.character(ID)
	b <- which(a == pattern[,1])
	pdf(paste0("Histogram modulus (", ID, ").pdf"), width = 11, height = 8) 
	hist(Data2[,4], xlim=c(0,1), breaks=seq(0,1,by=0.1), main = paste0("Exocytosis density vs normalized modulus (", ID, ") (n=", nrow(Data2), ")"), ylab="Density", xlab="Normalized modulus", col="firebrick", lwd=4)
	rect((1-pattern_radius), 0, 1, 1000, col= rgb(0.06,0.306,0.545,alpha=0.40), border = NA)
	dev.off()
	
	##plot Ripley's K function (single cell)
	if(Ripley==TRUE){
		lyso <- ppp(Data[,2]/radius,Data[,3]/radius,window=cellwindow)
		Kfunction <- envelope(lyso, Kest, nsim=100, nrank=1, transform = expression(. -pi*r^2), correction="best", savefuns=TRUE)
		pdf(paste0("Ripley's K function (", ID, ").pdf"), width = 11, height = 8)
		plot(Kfunction, xlab="Normalized radius (r)", ylab=expression(paste("K(r)-", pi, "r²")), main=paste0("Ripley's K function (", as.character(ID),")"), lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, legend=FALSE )
		dev.off()
		KfunctionList[[i]] <- Kfunction
	}

	
	#Put all exocytosis event in the same dataframe
	V1 <- as.vector(rep(ID, nrow(Data2)))
	Data2 <- cbind(Data2, V1)
	colnames(Data2) <- c("time", "Theta", "A", "A normalized", "ID")
	colnames(META) <- c("time", "Theta", "A", "A normalized", "ID")
	ifelse(nrow(META) < 2, META <- Data2, META <- rbind(META, Data2))
}

AveragedKfunction <- pool.envelope(KfunctionList[[1:length(KfunctionList)]])
pdf(paste0("Averaged Ripley's K function (N=", length(KfunctionList)," WT cells).pdf"), width = 11, height = 8)
plot(AveragedKfunction, xlab="Normalized radius (r)", ylab=expression(paste("K(r)-", pi, "r²")), main=paste0("Averaged Ripley's K function (N=", length(KfunctionList)," WT cells)"), lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, legend=FALSE )
dev.off()

##########################################################################################


##########################################################################################
## 5/ CSR test:


NND_MC <- vector(, N)
pvalue <- rep(-1, n)
Clustering <- 0
Dispersing <- 0


for(j in 1:n){
	A1 <- META[which(as.character(META[,5])==parameter[j,1]),]
	NND_Observed <- vector(, nrow(A1))
	if(nrow(A1) > 0) {
		A1 <- A1[,2:4]
		A1 <- A1[,-2]
		A2 <- A1
		colnames(A2) <- c("x", "y")
		A2[,1] <- cos(A1[,1]) * A1[,2]
		A2[,2] <- sin(A1[,1]) * A1[,2]
	
		for(i in 1:nrow(A1)){
			index <- 1:nrow(A1)
			index <- index[-i]
			NND_Observed[i] <- min(sqrt((A2[i,1] - A2[index,1])^2 + (A2[i,2] - A2[index,2])^2))
		}
		B <- mean(NND_Observed)

		#Monte-Carlo CSR simulation
		number_exo <- nrow(A2)
		for(k in 1:N) {	
			random_cell <- matrix(, ncol=2, nrow=number_exo)
	
			#Cell simulation
			for(i in 1:number_exo) {
				X <- 2
				Y <- 2
				while((X^2 + Y^2)>=1){
					Y <- runif(1, min = -1, max = 1)
					X <- runif(1, min = -1, max = 1)
				}
				random_cell[i,1] <- X
				random_cell[i,2] <- Y
		
			}
	
			NND_simulated <- vector(,number_exo)
			for(i in 1:number_exo){
				index <- 1:number_exo
				index <- index[-i]
				NND_simulated[i] <- min(sqrt((random_cell[i,1] - random_cell[index,1])^2 + (random_cell[i,2] - random_cell[index,2])^2))
			}	
	
			A <- mean(NND_simulated)
			NND_MC[k] <- A
	
		}
	
		NLess <- which(NND_MC < B)
		NLess <- length(NLess) + 1
		NGreater <- which(NND_MC > B)
		NGreater <- length(NGreater) + 1
	
		pvalueLess <- NLess/(length(NND_MC)+1)
		pvalueGreater <- NGreater/(length(NND_MC)+1)
		pvalue[j] <- min(pvalueLess,pvalueGreater)*2	#two-sided
		
		AverageNND_MC <- mean(NND_MC)
		ifelse(B>AverageNND_MC, Dispersing <- Dispersing +1, Clustering <- Clustering +1)
		
		pdf(paste0("NND distribution (", parameter[j,1], ").pdf"), width = 11, height = 8)
		d1 <- density(NND_MC, bw="SJ", kernel = "gaussian")
		plot(d1, main= paste0("NND distribution from Monte-Carlo simulation (p-value=", round(pvalue[j],3), ")"), ylab="Density", xlab="NND", col="blue", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
		abline(v = B, col="red", lwd=3, lty=2)
		dev.off()
	}
}

for(i in length(pvalue):1) {if(pvalue[i]<0) {pvalue <- pvalue[-i]}}

print(paste0("Clustering coefficient :", Clustering/length(pvalue)))
print(paste0("Dispersion coefficient :", Dispersing/length(pvalue)))

ks.test(pvalue,"punif",0,1)

pdf("P-value histogram of CSR test (Monte-Carlo).pdf", width = 11, height = 8) 
hist(pvalue, xlim=c(0,1), breaks=seq(0,1,by=0.1), main = "P-value histogram of CSR test (Monte-Carlo)", ylab="Count", xlab="P-value", col="cornflowerblue", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
mtext(paste0("Clustering coefficient :", round(Clustering/length(pvalue),3), "        Dispersion coefficient :", round(Dispersing/length(pvalue),3), "        p-valueKS =", round(ks.test(pvalue,"punif",0,1)$p,3)), cex=1.3)
dev.off()
##########################################################################################


##########################################################################################
## 6/ Random sampling from original datasets:

if(RandomSampling == TRUE ) {
	META2 <- as.data.frame(matrix(,nrow=0,ncol=5))
	colnames(META2) <- colnames(META)
	Effectif <- vector(,length=nrow(pattern))
	for(i in 1:nrow(pattern)){
		A1 <- which(pattern[i,1] == META[,5])
		Effectif[i] <- length(A1)
	}

	for(i in 1:nrow(pattern)){
		A1 <- which(pattern[i,1] == META[,5])
		A2 <- sample(A1, min(Effectif), replace=FALSE)
		META2 <- rbind(META2, META[A2,])
	}
}

if(RandomSampling == FALSE){
	META2 <- META
}

##########################################################################################


##########################################################################################
## 7/ Modulus distribution:

Envellope_Modulus <- matrix(0, nrow=(Echan+1), ncol=N_sim)
Envellope_FModulus <- matrix(0, nrow=(Echan+1), ncol=N_sim)

for(j in 1:N_sim) {
	##Monte-Carlo CSR simulation
	NbSimulatedCell <- length(Effectif)*min(Effectif)
	random_cell <- as.data.frame(matrix(, ncol=2, nrow=NbSimulatedCell))
	colnames(random_cell) <- c("modulus", "Angle")
	for(i in 1:NbSimulatedCell) {
		X <- 2
		Y <- 2
		while((X^2 + Y^2)>=1){
			Y <- runif(1, min = -1, max = 1)
			X <- runif(1, min = -1, max = 1)
		}
		random_cell[i,1] <- sqrt(X^2 + Y^2)
		if(Y >= 0) {random_cell[i,2] <- acos(X/random_cell[i,1]) * 180/pi}
		if(Y < 0) {random_cell[i,2] <- acos(-X/random_cell[i,1]) * 180/pi +180}
	}

	d1 <- dbckden(seq(0,1,1/Echan),random_cell[,1], lambda = 0.04, kernel = "gaussian", bcmethod = "beta1", proper = TRUE, nn = "jf96", offset = NULL, xmax = 1, log = FALSE)
	Envellope_Modulus[,j] <- d1
	
	p <- ecdf(random_cell[,1])
	for(i in 0:(1*Echan)){
		Envellope_FModulus[i,j] <- p(i/Echan)
	}
}

#Monte Carlo Envellope (1 and 99 percentile)
Envellope_Modulus2 <- matrix(0, ncol=3, nrow=(Echan+1))
colnames(Envellope_Modulus2) <- c("Modulus", "Min value", "Max value")
Envellope_Modulus2[,1] <- seq(0,1, (1/Echan))
for(i in 1:(Echan+1)) {Envellope_Modulus2[i,2] <- quantile(Envellope_Modulus[i,], prob = 1-alpha)}
for(i in 1:(Echan+1)) {Envellope_Modulus2[i,3] <- quantile(Envellope_Modulus[i,], prob = alpha)}

PatternRadius <- mean(pattern[,4])
SEMPatternRadius <- sd(pattern[,4])/sqrt(nrow(pattern))

#Plot
pdf("Modulus distribution.pdf", width = 11, height = 8) 
p <- dbckden(seq(0,1,0.001),META2[,4], lambda = 0.04, bcmethod = "beta1", proper = TRUE, nn = "jf96", offset = NULL, xmax = 1, log = FALSE)
plot(y=p, x=seq(from = 0, to =1, by =0.001), type="l", main= paste0("Modulus distribution (n=", length(Effectif), "x", min(Effectif), ")"), ylab="Density", xlab="Normalized modulus", col="red", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, axes=FALSE)
rect((1-PatternRadius), 0, 1, 1000, col= adjustcolor("darkgray",alpha=0.40), border = NA)
axis(side=2, at=graduation4, labels=vec.expr2, cex.axis = 1.2, lwd=4)
axis(side=1, at=graduation2, labels=vec.expr1, cex.axis = 1.2, lwd=4)
par(new=TRUE)
curve(2*x, col="blue", lwd=4, lty=3, xaxt='n', yaxt='n', ann=FALSE,  ylim=c(0,2))
par(new=FALSE)
errorbar((1-PatternRadius), 0.5, (1-PatternRadius-SEMPatternRadius), (1-PatternRadius+SEMPatternRadius), 0.2)

#Error band at 99%
xshade1 <- c(Envellope_Modulus2[,1],rev(Envellope_Modulus2[,1]))
yshade1 <- c(Envellope_Modulus2[,2], rev(Envellope_Modulus2[,3]))
polygon(xshade1,yshade1, border = NA,col=adjustcolor("dodgerblue2", 0.4))

estimates <- matrix(NA,nrow=(evalpoints+1),ncol=NBoot)
for (b in 1:NBoot){
  xstar <- sample(META2[,4],replace=TRUE)
  dstar <- dbckden(seq(0,1,1/evalpoints),xstar, lambda = 0.04, kernel = "gaussian", bcmethod = "beta1", proper = TRUE, nn = "jf96", offset = NULL, xmax = 1, log = FALSE)
  estimates[,b] <- dstar
}
ConfidenceBands <- apply(estimates, 1, quantile, probs = c(1-alpha, alpha))
xax <- seq(0,1,length.out = evalpoints+1)
xshade <- c(xax,rev(xax))
yshade <- c(ConfidenceBands[2,],rev(ConfidenceBands[1,]))
polygon(xshade,yshade, border = NA,col=adjustcolor("firebrick", 0.4))
dev.off()

##########################################################################################
## 8/ CDF:

CDF_Modulus <- function(x){
	if(x>= 0 && x<=1) {return(x^2)}
	if(x<0) {return(0)}
	if(x>1) {return(1)}
}
CDF_Modulus <- Vectorize(CDF_Modulus)

CDF <- ecdf(META2[,4])

pdf("Modulus empirical distribution.pdf", width = 11, height = 8) 
plot(CDF, main= paste0("Exocytosis modulus empirical distribution (n=", length(Effectif), "x", min(Effectif), ")"), ylab="F(modulus)", xlab="Normalized modulus", col="red", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, xlim=c(-0.05,1.05), axes=FALSE)
rect((1-PatternRadius), 0, 1, 1000, col= adjustcolor("darkgray",alpha=0.40), border = NA)
axis(side=2, at=graduation4, labels=vec.expr2, cex.axis = 1.2, lwd=4)
axis(side=1, at=graduation2, labels=vec.expr1, cex.axis = 1.2, lwd=4)
par(new=TRUE)
curve(CDF_Modulus, col="blue", lwd=4, lty=3, xaxt='n', yaxt='n', ann=FALSE,  ylim=c(0,1), xlim=c(-0.05,1.05))
par(new=FALSE)
errorbar((1-PatternRadius), 0.2, (1-PatternRadius-SEMPatternRadius), (1-PatternRadius+SEMPatternRadius), 0.1)

Envellope_FModulus2 <- matrix(0, ncol=3, nrow=(Echan+1))
colnames(Envellope_FModulus2) <- c("F(Modulus)", "Min value", "Max value")
Envellope_FModulus2[,1] <- seq(0,1, (1/Echan))
for(i in 1:(Echan+1)) {Envellope_FModulus2[i,2] <- quantile(Envellope_FModulus[i,], prob = 1-alpha)}
for(i in 1:(Echan+1)) {Envellope_FModulus2[i,3] <- quantile(Envellope_FModulus[i,], prob = alpha)}
xshade1 <- c(Envellope_FModulus2[,1],rev(Envellope_FModulus2[,1]))
yshade1 <- c(Envellope_FModulus2[,2], rev(Envellope_FModulus2[,3]))
polygon(xshade1,yshade1, border = NA,col=adjustcolor("dodgerblue2", 0.4))

estimates <- matrix(NA,nrow=(evalpoints+1),ncol=NBoot)
for (b in 1:NBoot){
  xstar <- sample(META2[,4],replace=TRUE)
  dstar <- ecdf(xstar)
  for(j in 1:(evalpoints+1)){
	estimates[j,b] <- dstar(j/(evalpoints+1))
  }
}
ConfidenceBands <- apply(estimates, 1, quantile, probs = c(1-alpha, alpha))
xax <- seq(0,1,length.out = evalpoints+1)
xshade <- c(xax,rev(xax))
yshade <- c(ConfidenceBands[2,],rev(ConfidenceBands[1,]))
polygon(xshade,yshade, border = NA,col=adjustcolor("firebrick", 0.4))
dev.off()