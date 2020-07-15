# simulate a dataset
###############INPUTS and Outputs#########################

testRun <- FALSE

if (testRun == TRUE) {
	mainPath <- "/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/"
	outPutPathPlot <- paste0(mainPath, "ZS/SimCreation.pdf")
	outPutPathSimCompDF <- paste0(mainPath, "ZS/Sim_completeDF.rda")
	outPutPathPropDF <- paste0(mainPath, "ZS/SimPropDF.rda")
}
if (testRun == FALSE) {
	outPutPathPlot = snakemake@output[["plotCreat"]]
	outPutPathSimCompDF = snakemake@output[["SimCompVar"]]
	outPutPathPropDF = snakemake@output[["propDFVar"]]
}


###############LIBRARIES#########################

library(stringr)
library(brms)

#################CODE############################




geneList <- paste0("GeneID_",str_pad(seq(1,20000), 5, pad = "0"))
groupList <- rep(0,length(geneList))
groupList[0:331] <- 1

ListOfMutations <- c()
for (n in letters[1:15]) {
	for (m in letters[1:15]) {
		ListOfMutations <- c(ListOfMutations, paste0(n,".",m))
	}
}


	
	
BGProp <- runif(length(ListOfMutations), min=0,max = 1)
GOProp <- c()
for (n in BGProp) {
	
	
	check = TRUE
	while (check == TRUE) {
		x <- rnorm(1,mean=n, sd= 0.2)
		
		if (x < 1 & x > 0) { check = FALSE	}
		
	}
	GOProp <- c(GOProp, x)
	
	
}
#rm(x)

propDF <- data.frame(mut=ListOfMutations,BGProp = BGProp, GOProp =GOProp)
rownames(propDF) <- propDF$mut

pdf(outPutPathPlot)
hist((BGProp) , main= "Draws for 225 Mutations", xlab= "Probability of Occurrence in BG")
hist((GOProp), col= "red", main= "Draws for 225 Mutations", xlab= "Probability of Occurrence in Gene-Subset")


plot(density(logit_scaled(BGProp) - logit_scaled(GOProp)), col= "blue", main= "Distribution of Enrichment and Depletion" , xlab= "logOdd-Ratio between BG and Gene-Subset")
dev.off()





#############Generate

bigDF <- as.data.frame(matrix(NA, nrow = length(geneList), ncol = 2+length(ListOfMutations)))
rownames(bigDF) <- geneList
colnames(bigDF) <- c("Gene","Group", ListOfMutations)

bigDF$Gene <- geneList
bigDF$Group <- groupList


for (mut in ListOfMutations) {
	valueVector <- c(rbinom(334, 1, propDF[mut,"GOProp"]),rbinom(length(rownames(bigDF))-334, 1, propDF[mut,"BGProp"]))
	
	bigDF[,mut] <- valueVector
}


###########Prep for read in
simCompletDFA <- as.data.frame(matrix(NA, nrow = length(ListOfMutations), ncol = 4))
colnames(simCompletDFA) <- c("mut"  ,"hits" ,"trails", "GO")
simCompletDFA$mut <- ListOfMutations
simCompletDFB <- simCompletDFA

simCompletDFA$GO <- 0
simCompletDFB$GO <- 1



for (mut in ListOfMutations) {
	x <- (bigDF[bigDF$Group == 0, mut])
	y <- (bigDF[bigDF$Group == 1, mut])
	
	simCompletDFA[simCompletDFA$mut == mut, "hits"] <- sum(x)
	simCompletDFA[simCompletDFA$mut == mut, "trails"] <- length(x)
	
	simCompletDFB[simCompletDFB$mut == mut, "hits"] <- sum(y)
	simCompletDFB[simCompletDFB$mut == mut, "trails"] <- length(y)
}

simCompletDF <- rbind(simCompletDFA,simCompletDFB)

save(propDF, file=outPutPathPropDF)
save(simCompletDF, file=outPutPathSimCompDF)