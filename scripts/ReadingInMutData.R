#Reading in the chromHMM data (state transitions mapped to the genes)
###############INPUTS and Outputs#########################
#rawInputDataPath <- "/media/ubuntuDisk/home/knalli/Dropbox/Phd/PhD-Work_Drop/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment/data/Jan2020gcsc3.mutations.RData"
#threeGOPath      <- "/media/ubuntuDisk/home/knalli/Dropbox/Phd/PhD-Work_Drop/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment/data/ListOfGenes/HKGList.rda"

x = load(snakemake@input[["mutHitDF"]])
#x <- load(rawInputDataPath)
raw <- get(x) 
rm(x)

#x <- load(threeGOPath)
x = load(snakemake@input[["GeneSelection"]])
GOList <- unique(get(x) )
rm(x)

#################CODE############################
BGList <- unique(raw$gene)
regionName <- c("enhancer", "tss"  ,    "gene")[3]



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
completeDF <- data.frame(mut=compList,PIndex= PIntigrate, HIndex= HListIntigrate)




for (n in c("inBG","inGO")) {
	geneDF <- NA
	if (n == "inBG") {geneDF <- raw[raw$gene %in% BGList & raw$region == regionName, ]			}
	if (n == "inGO") {geneDF <- raw[raw$gene %in% GOList & raw$region == regionName, ]			}
	
	#geneDF <- dupFunc(geneDF) #islandcount includes
	#return(geneDF)
	temp <- as.data.frame(table(geneDF[ , c("state")]))
	colnames(temp) <- c("mut", paste0("hits"))
	
	#temp$trails <- length(unique(geneDF$gene)	)
	

	if (n == "inBG") {tempBG  <- merge(completeDF, temp, by = "mut", all.x = TRUE)
		tempBG$GO <- 0
		tempBG$trails <- length(unique(geneDF$gene)	)
	}
	if (n == "inGO") {tempGO <-  merge(completeDF, temp, by = "mut", all.x = TRUE)	
		tempGO$GO <- 1
		tempGO$trails <- length(unique(geneDF$gene)	)
	}	
	
	
	
}

completeDF <- rbind(tempBG,tempGO)
completeDF[is.na(completeDF)] <- 0

completeDF <- completeDF[,c("mut", "hits", "trails", "GO")]

save(completeDF, file = snakemake@output[["compDFName"]])




