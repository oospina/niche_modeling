## Set output folder
outputfp <- "C:/Users/OSCAREOSPINATOBON/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_MaxEntModels_Feb042020/maxEntModel_Axera_Feb042020_Plots"

## Set path of "plots" folder within the MaxEnt run
inputfp <- "C:/Users/OSCAREOSPINATOBON/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_MaxEntModels_Feb042020/maxEntModel_Axera_Feb042020/plots"

## Set species name as indicated in the html MaxEnt output
species <- "xera"

## Set number of MaxEnt replicates
nrep <- 30

## Set name of variable to be plotted as inputted to MaxEnt without the ".asc" suffix
variable <- "bio_aligned_14"

## Set axis name to be displayed on X axis
varname <- "Precipitation of Driest Month"

## Set type of variable : "CONT" OR "CAT"
vartype <- "CONT"

####################
if(vartype == "CONT" || vartype == "CAT") {
	vartype <- vartype
}else {
	print("ERROR: Select \"CONT\" or \"CAT\".", quote=FALSE)
}

library("Hmisc")
setwd(outputfp)

xvalues <- read.csv(paste(inputfp, "/", species, "_0_", variable, "_only.dat", sep=""))[2] ## Takes values for X axis
df <- data.frame(xvalues)

for(i in 0:(nrep-1)) {
	yvalues_tmp <- read.csv(paste(inputfp, "/", species, "_", i, "_", variable, "_only.dat", sep=""))[3]
	df <- cbind(df, yvalues_tmp)
	colnames(df)[i+2] <- paste("y", i, sep="")	
}

varmeans <- rowMeans(df[, 2:(nrep+1)])
varminus <- varmeans - (apply(df[, 2:nrep+1], 1, sd))
varplus <- varmeans + (apply(df[, 2:nrep+1], 1, sd))
df <- cbind(df, varmeans, varminus, varplus)

if(vartype == "CONT") {
	print(min(df[(df$varmeans>=0.5), ]$x))
	print(max(df[(df$varmeans>=0.5), ]$x))
}

if(vartype == "CONT") {

	pdf(paste(outputfp, "/", species, "_", variable, ".pdf", sep=""))
	plot(df$x, df$varmeans, ylim=c(0,1), type="n", ylab="Probabilty of ocurrence", xlab=varname, xlim=c(min(df$x),max(df$x)))
	lines(df$x, df$varmeans, lty=1, lwd=3)
	lines(df$x, df$varminus, lty=1, col="black")
	lines(df$x, df$varplus, lty=1, col="black")
	abline(h=0.5, lty=2, col="gray")
	dev.off()
}

if(vartype == "CAT") {

	df_sub <- df[order(df$varmeans, decreasing=T), ]
	df_sub <- df_sub[1:20, ]
	
	pdf(paste(outputfp, "/", species, "_", variable, ".pdf", sep=""))
	par(mar=c(6, 4.1, 4.1, 2.1))
	par(mgp=c(4.5,1,0))
	barplot(df_sub$varmeans, ylim=c(0,1), ylab="Probabilty of ocurrence", xlab=varname, names=df_sub$x, col="white", las=2)
#	bp <- barplot(df_sub$varmeans, plot=F)
#	errbar(bp, df_sub$varmeans, df_sub$varplus, df_sub$varminus, add=T)
	abline(h=0.5, lty=2, col="gray")
	dev.off()

}
rm(list=ls())
