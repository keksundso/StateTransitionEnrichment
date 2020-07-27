#get posterior prediction of brms model, and plotting it
###############INPUTS and Outputs#########################
testRun <- FALSE

if (testRun == TRUE) {
  parentPath <- "/media/knalli/07eee2b2-1320-41f3-b7bb-336496e48287/home/knalli/Phd/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment"
	x <- load(file.path(parentPath,"ZS/HKGList_completeDF.rda"))
	y <- load(file.path(parentPath,"ZS/HKGList_BinomialModle.rda"))

	outPutPath <- (file.path(parentPath,"/testPlays/prettyPPC_20200429.tiff"))

}
if (testRun == FALSE) {
	x = load(snakemake@input[["compDFName"]])
	y = load(snakemake@input[["brmsModleVar"]])
	outPutPath = snakemake@output[["PPCVar"]]
}

completeDF <- get(x)
myModel <- get(y) 
rm(x)

###############LIBRARIES#########################

library(brms)
myMod <- function(valueVector){
	d <- density(valueVector)
	df <- data.frame(d$x,d$y)
	
	value <- df[df$d.y == max(df$d.y), "d.x" ]
	
	return(value)
}

#################CODE############################

Ueberschrift <- "PPC of the HKG-Selection Model"
if (testRun == FALSE) {pdf(outPutPath, width=1500, height=1300)
	Ueberschrift <- paste("PPC:",basename(snakemake@input[["brmsModleVar"]]))
	}


completeDFBG <- completeDF[completeDF$GO == 0 ,]
completeDFGO <- completeDF[completeDF$GO == 1 ,]

listOFmut <- completeDFBG[order(completeDFBG$hits),"mut"]
#listOFmut <- completeDFBG$mut

modleNameList <- c("Binomial", "Beta Binomial")
outerCounter <- 0


postPredict <- posterior_predict(myModel)
colnames(postPredict) <- c(paste(completeDFBG$mut,"BG",sep="_"),paste(completeDFGO$mut,"GO",sep="_"))



plot(c(-1000,-225),c(-1,-500),xlim=c(0,225),ylim=c(0,max(completeDFBG$trails))
	 , xlab= "transition index", ylab= "occurrence", main = Ueberschrift)


counter <- 0
for (mut in listOFmut) {
	counter <- counter +1
	
	y <- postPredict[,paste0(mut,"_BG")]
	tmp <-  quantile(y, probs = c(0.11,0.89))
	#print(paste0(counter,completeDFBG[completeDFBG$mut == mut,"hits" ]))
	points(counter,completeDFBG[completeDFBG$mut == mut,"hits" ],pch=19,cex=0.5, col="black")
	lines(c(counter-0.75,counter+0.75), c(myMod(y),myMod(y)), col= "gray",lwd= 1.5)
	lines(c(counter,counter), c(tmp[1]  ,tmp[2]), col= "gray")
}


par(new = TRUE)
plot(c(-1000,-225),c(-1,-500),xlim=c(0,225),ylim=c(-20,max(completeDFGO$trails)),xaxt="n", yaxt="n"
	 , xlab= "", ylab= "", main = ""
)
ytick<-seq(0, max(completeDFGO$trails), by=75)


axis(side=4, at=ytick, labels = FALSE)
text(par("usr")[2], ytick,  col="#7CAE00",
	 labels = ytick, srt = 0, pos = 4, xpd = TRUE)

counter <- 0
for (mut in listOFmut) {
	x <- postPredict[,paste0(mut,"_GO")]
	counter <- counter +1
	factorValue <- 1
	
	tmp <-  quantile(x, probs = c(0.11,0.89))
	
	points(counter,completeDFGO[completeDFGO$mut == mut,"hits" ]*factorValue,pch=19, col="#7CAE00",cex=0.5)
	lines(c(counter-0.75,counter+0.75), c(myMod(x)*factorValue      ,myMod(x)*factorValue), col= "gray",lwd= 1.5)
	lines(c(counter,counter),         c(tmp[1]  ,tmp[2]), col= "gray")
}
legend("topleft",legend = c("In the Background","In the HKG-List","Posterior" ),col=c("black","#7CAE00","gray"), pch=19)
# Legend text does not adapt to the input data (will be fixed if moved to ggplot)

if (testRun == FALSE) {dev.off()}




