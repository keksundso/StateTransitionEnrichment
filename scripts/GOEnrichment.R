#Perform an Enrichment analysis (not just GO) on the subset of genes which are in the group of interest and show the transition
#the background for the analysis are all genes in the group of interest.
###############INPUTS and Outputs#########################
testRun <- FALSE

if (testRun == TRUE) {
	mainPath <- "/media/ubuntuDisk/home/knalli/Dropbox/Phd/PhD-Work_Drop/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment/"
	
	x <- paste0(mainPath,"data/Jan2020gcsc3.mutations.RData")
	y <- paste0(mainPath,"data/ListOfGenes/GOthreeComb.rda")
	z <- paste0(mainPath,"ZS/GOthreeComb_BinomialModle.rda")
	
	outPutPath <- paste0(mainPath,"testPlays/GOAnalysis20200512/GOEnrichment.pdf")
}
if (testRun == FALSE) {
	x = snakemake@input[["mutHitDF"]]
	y = snakemake@input[["GOVector"]]
	z = snakemake@input[["brmsModleVar"]]  
	outPutPath = snakemake@output[["pdfEnrichment"]]
	
	print(x)
	print(y)
	print(outPutPath)
}

inputDataTabel <- get(load(x))
geneOfIntrest <- get(load(y) )
brmsModle <- get(load(z) )
rm(x,y,z)

###############LIBRARIES#########################

#install.packages("devtools")
#library("devtools")
#devtools::install_github("yingstat/goview", dependencies = FALSE)
library(goview)
library(gProfileR)
library(ggplot2)
library(pdftools)
library(brms)


#################CODE############################


myMod <- function(valueVector){
  #Return most likely value of density distribution
	d <- density(valueVector)
	df <- data.frame(d$x,d$y)
	
	value <- df[df$d.y == max(df$d.y), "d.x" ]
	
	return(value)
}

extractingDataFromFunction <- function(brmsModel){
	listOfMut <- c()
	listOfHCCol <- c()
	listOfPCType <- c()
	
	PCseq <- seq(1,15)
	HCseq <- seq(1,15)
	
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
		zeroProp <- pnorm(0,mean=mean(y),sd(y))
		
		#print(paste(n, length(y[y<=0])))
		
		completeDF[realMut,"modalWert"] <- modalWert
		completeDF[realMut,"zeroProp"] <- zeroProp
		
		if (modalWert <= 0) {ZeroCount <- length(y[y>=0]) }
		if (modalWert >= 0) {ZeroCount <- length(y[y<=0]) }
		
		completeDF[realMut,"countOfZero"] <- (ZeroCount/ length(y))   +(1/length(y))
		
	}
	
	normProp <- c()
	for(row in rownames(completeDF)){
		if(completeDF[row,"zeroProp"] < 0.5){normProp <- c(normProp,  completeDF[row,"zeroProp"])}
		if(completeDF[row,"zeroProp"] > 0.5){normProp <- c(normProp, 1-   completeDF[row,"zeroProp"])}
		if(completeDF[row,"zeroProp"] == 0){print(row)}
	}
	completeDF$normProp <- normProp
	
	return(list(completeDF, length(y)))
}

#extract genes which have a mutation in a given set

getDataFrameForGo <- function(outputPath ,m, geneListOfIntrest, mutOfIntrest){

	
	dir.create(outputPath,showWarnings = FALSE)
	
	compList <- c()
	PIntigrate <- c()
	HListIntigrate <- c()
	for (P in seq(1:15)) {
		for (H in seq(1:15)) {
			compList <- c(compList, paste0(P,".",H))
			PIntigrate <- c(PIntigrate, P)
			HListIntigrate <- c(HListIntigrate, H)
		}
	}
	completeDF <- data.frame(genes = geneListOfIntrest)
	for (mut in compList) {
		completeDF <- cbind(completeDF, mut= rep(0,n=length(geneListOfIntrest)))
		colnames(completeDF) <- c(head( colnames(completeDF), -1),mut )
	}
	rownames(completeDF) <- completeDF$genes
	#completeDF$genes <- NULL
	completeDF$backGround <- 1
	
	mForGenes <- m[m$gene %in% geneListOfIntrest & m$region == c("gene"), ]
	
	for (line in rownames(mForGenes)) {
		completeDF[as.character(mForGenes[line, "gene"]) , as.character(mForGenes[line, "state"]) ] <- 1
	}
	
	#	return(completeDF)
	counter <- 0
	for (colName in colnames(completeDF)) {
		counter <- counter + 1
		#if (counter%%20 == 0 ) {print(counter)}
		colNameA <- gsub('\\.', '-', colName)
		tempGeneList <- completeDF[completeDF[,colName]==1, "genes"]
		tempGeneListA <- c()
		for (n  in tempGeneList) {
			tempGeneListA <- c(tempGeneListA, n)
		}
		
		
		#if (length(tempGeneList) > 0 ) {
		if (colName %in% mutOfIntrest){
			x <- gprofiler(tempGeneListA, organism = "mmusculus",custom_bg = as.vector(completeDF$genes ))
			#browser()
			save(tempGeneListA, file=paste0(dirname(outPutPath),"/tempGeneListA.Rda"))
			save(completeDF, file=paste0(dirname(outPutPath),"/completeDF.Rda"))
			if (length(x$term.size) < 1) {	colNameA <- paste0("n_",colNameA)	}
			#pdf(paste0(outputPath,"/",colNameA,".pdf"))
			y <- goview(x)
			y$labels$title <- colNameA
			save(x, file = paste0(outputPath,"/",colNameA,"goviewInput.Rda") )
			save(y, file = paste0(outputPath,"/",colNameA,".ggplot") )
			ggsave(filename =paste0(outputPath,"/",colNameA,".pdf"), y, width = 30, height = 15,units = "cm")
			write(tempGeneListA, file = paste0(outputPath,"/",colNameA,".txt"), sep = "\n")
			write.csv2(x,             file = paste0(outputPath,"/",colNameA,"_GOEnr.csv"))
			
		}
	}
	
	return(y)
}

dataFromModle <- extractingDataFromFunction(brmsModle)
mutOfIntrest <- unique(head(as.character(dataFromModle[[1]]$mut[order(dataFromModle[[1]]$normProp)]), n= 25))

tempListDF<- getDataFrameForGo(dirname(outPutPath),inputDataTabel, geneOfIntrest, mutOfIntrest)


allPDFs <- list.files(dirname(outPutPath), pattern = c(".pdf"))
unwantedPDFs <- c(list.files(dirname(outPutPath), pattern = c("n_")),"GOEnrichment.pdf")
wantedPDFs <- setdiff(allPDFs,unwantedPDFs)
wantedPDFs <- file.path(dirname(outPutPath),setdiff(allPDFs,unwantedPDFs))

pdf_combine(wantedPDFs ,outPutPath)
