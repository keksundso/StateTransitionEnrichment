#creating the transition enrichment plots ("probability of direction" against log(odds ratio))
###############INPUTS and Outputs#########################

testRun <- FALSE

if (testRun == TRUE) {
  parentPath <- "/media/ubuntuDisk/home/knalli/Dropbox/Phd/PhD-Work_Drop/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment"
	ZSPath <- file.path(parentPath,"ZS/")
	x <- paste0(ZSPath, "GOthreeComb_completeDF.rda")
	y <-paste0(ZSPath, "GOthreeComb_BinomialModle.rda")
	
	outPutPath <- file.path(parentPath,"testPlays/testOutputPlot.pdf")
}

if (testRun == FALSE) {
	x = snakemake@input[["compDFName"]]
	y = snakemake@input[["brmsModleVar"]]
	outPutPath = snakemake@output[["plotVar"]]
	print(x)
	print(y)
	print(outPutPath)
}

completeDF <- get(load(x))
meinModel <- get(load(y) )
rm(x)

###############LIBRARIES#########################

library(brms)
normPropPlot <- TRUE # Probability that a normal distribution would cross zero and not the empirical posterior


myMod <- function(valueVector){
	d <- density(valueVector)
	df <- data.frame(d$x,d$y)
	
	value <- df[df$d.y == max(df$d.y), "d.x" ]
	
	return(value)
}

#################CODE############################

# Plotting the Modle
#- x-Axis is the mean of the posterior Slope for each transition
#- y-Axis is the percentage of the Posterior Distribution which is above Zero (or below in the case of a depletion)

##plotParameter
listOfActivStates <- c(4,5,7,8,9,10,11,12)
listOfBiStates <- c(3,13,15)
listOfRepStates <- c(1,2,14)
listOfEmtpyStates <- c(6)

PCType <- c(15,18,19,7)
HCCol <- c("green","yellow","red","grey")

plotAES <- data.frame(state = c(listOfActivStates,listOfBiStates,listOfRepStates,listOfEmtpyStates),
					  PCType = c(rep(15, length(listOfActivStates)),rep(18, length(listOfBiStates)),rep(19, length(listOfRepStates)),rep(7, length(listOfEmtpyStates))  ) ,
					  HCCol = c(rep("green", length(listOfActivStates)),rep("yellow", length(listOfBiStates)),rep("red", length(listOfRepStates)),rep("grey", length(listOfEmtpyStates))  ) 
)

rownames(plotAES) <- plotAES$state

listOfMut <- c()
listOfHCCol <- c()
listOfPCType <- c()

PCseq <- seq(1,15)
HCseq <- seq(1,15)
if (testRun == FALSE) {
	if (grepl("Sim_completeDF.rda", snakemake@input[["compDFName"]]) ) {
		PCseq <- letters[1:15]
		HCseq <- letters[1:15]
		print("Sim Complete")
	}
	if (!(grepl("Sim_completeDF.rda", snakemake@input[["compDFName"]])) ) {
		PCseq <- seq(1,15)
		HCseq <- seq(1,15)
		print("CompleteDF no sim")
	}
}

for (PC in PCseq) {
	for (HC in HCseq) {
		listOfMut    <- c(listOfMut, paste0(PC,".",HC))
		listOfPCType <- c(listOfPCType, plotAES[as.character(PC), "PCType"])
		listOfHCCol  <- c(listOfHCCol, as.character(plotAES[as.character(HC), "HCCol"]))
	}
}
completeDFAES <- data.frame(mut=listOfMut, PCType = listOfPCType, HCCol = listOfHCCol)	

rm(listOfActivStates,listOfBiStates, listOfRepStates,listOfEmtpyStates, PCType, HCCol, plotAES, listOfMut,listOfHCCol, listOfPCType )




#load(file=paste0(pathParameter,"ZS/brmsModle_",RunName, ".Rda"))


modlePlotFunction <- function(brmsModel){
  # Function to create the data table which contains the values for the plot
	x <- (posterior_samples(brmsModel))
	
	#y <- x$`r_mut[11.15,GO]`
	
	#plot(c(-100,100),c(100,199), xlim=c(-3,3),ylim = c(0,2.5), xlab = "Effect Size", ylab = "Density")
	
	compList <- c()
	PIntigrate <- c()
	HListIntigrate <- c()
	for (P in PCseq) {
		for (H in HCseq) {
			compList <- c(compList, paste0(P,".",H))
			PIntigrate <- c(PIntigrate, P)
			HListIntigrate <- c(HListIntigrate, H)
		}
	}
	
	completeDF <- data.frame(mut=compList,PIndex= PIntigrate, HIndex= HListIntigrate)
	completeDF$modalWert <- NA
	completeDF$zeroProp <- NA
	completeDF$countOfZero <- NA
	rownames(completeDF) <- completeDF$mut
	
	for(n in colnames(x)[grepl("GO",colnames(x))&grepl("r_mut\\[",colnames(x))]){
		y <- x[, n]
		
		#lines(density(y))
		#lines(density(rnorm(1000000,mean=mean(y),sd(y))), col="red")
		realMut <- NA
		for(mut in unique(completeDF$mut)){  if(grepl(mut,n)){realMut <- mut} }
		#print(paste(realMut, n))
		
		modalWert <- myMod(y)
		#print(modalWert)
		zeroProp <- pnorm(0,mean=mean(y),sd(y)) #probability that a normel distro would cross zero
		
		#print(paste(n, length(y[y<=0])))
		
		completeDF[realMut,"modalWert"] <- modalWert
		completeDF[realMut,"zeroProp"] <- zeroProp
		
		if (modalWert <= 0) {ZeroCount <- length(y[y>=0]) }
		if (modalWert >= 0) {ZeroCount <- length(y[y<=0]) }
		
		completeDF[realMut,"countOfZero"] <- (ZeroCount/ length(y))   +(1/length(y))
		
	}
	
	normProp <- c() # probability that a normal distro would cross zero, but only beloew
	for(row in rownames(completeDF)){
		if(completeDF[row,"zeroProp"] < 0.5){normProp <- c(normProp,  completeDF[row,"zeroProp"])}
		if(completeDF[row,"zeroProp"] > 0.5){normProp <- c(normProp, 1-   completeDF[row,"zeroProp"])}
		if(completeDF[row,"zeroProp"] == 0){print(row)}
	}
	completeDF$normProp <- normProp
	
	return(list(completeDF, length(y)))
}


if (testRun == FALSE) {pdf(outPutPath)}
par(pty="s")

Uberschrift <- paste0(basename(y))
u <- modlePlotFunction(meinModel)
x <- u[[1]]
minValue <- 1/u[[2]]
rm(u)

x <- merge(x, completeDFAES, by="mut")
rownames(x) <- x$mut



if (normPropPlot == TRUE) {
	plot(x$modalWert,(x$normProp), log = "y"
		 , pch= as.numeric(x$PCType) , col = as.character( x$HCCol)
		 , xlab = "log(Odds-Ratio)", ylab= "Probability of no or opposite Effect (norm)", main=Uberschrift
		 )
	# very unelegant way to draw the line (will be moved to ggplot)
	lines(c(-0,0)  ,c(1,0.0000000000000000000000000000000000000000000000000000000000000001), lty= 3, col="gray")
	
	
	text(x$modalWert,(x$normProp), labels = rownames(x),cex = 0.2)
}

if (normPropPlot == FALSE) {
	plot(x$modalWert,(x$countOfZero), log = "y"
		 , pch= as.numeric(x$PCType) , col = as.character( x$HCCol)
		 , xlab = "log(Odds-Ratio)", ylab= "Probability of no or opposite Effect (empirical)", main=Uberschrift)
	
	lines(c(-110,100)  ,c(minValue,minValue), lty= 3, col="gray")
	lines(c(-0,0)  ,c(1,0.0000000000000000000000000000000000000000000000000000000000000001), lty= 3, col="gray")
	
	
	text(x$modalWert,(x$countOfZero), labels = rownames(x),cex = 0.2)
}

if (TRUE) {
	PCType <- c(15,18,19,7)
	HCCol <- c("green","yellow","red","grey")
	
	plot(c(-5,5),c(-5,5), col= "white", xlim = c(-5,5), ylim = c(-3,8)
		 ,xaxt='n',yaxt='n', xlab="", ylab = "",bty="n")
	points( c(-1,-1,-1,-1),c(4,3,2,1), pch= PCType)
	points( c(1,1,1,1),c(4,3,2,1), pch = 3, col= HCCol)
	
	text(-1,5, labels = "PC")
	text(1,5, labels = "HC")
	
	text(0,4, labels = "Activ")
	text(0,3, labels = "Bi")
	text(0,2, labels = "Rep")
	text(0,1, labels = "Empty")
}

if (testRun == FALSE) {dev.off()}


