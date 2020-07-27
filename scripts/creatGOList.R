#takes the names of GO-Terms (specified in snakemake (GOListName="data/GOTerms/{GOList}.go"))
# and outputs all genes associated to it

###############INPUTS and Outputs#########################
testRun <- FALSE

if (testRun == TRUE ) {
	inputPath <- "/Volumes/entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/data/GOTerms/GOgliogenesis.go"
	outputPath <- "/Volumes/entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/data/ListOfGenes/GOgliogenesis.rda"
}

if (testRun == FALSE ) {
	inputPath <- snakemake@input[["GOListName"]]
	outputPath <- snakemake@output[["GOVector"]]
}
###############LIBRARIES#########################
#integrated in the code section 

#################CODE############################
GOTERMListFunction <- function(GOID, dateiName="notDefined"){
	library(biomaRt)
	ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl"
					   #,host="apr2019.archive.ensembl.org"
	)
	if (dateiName == "notDefined") {
		GoTerm <- data.frame(getBM(attributes=c("name_1006","go_id"),
								   filters=c('go'),
								   values = GOID ,
								   mart=ensembl))
		
		dateiName <- GoTerm[GoTerm$go_id ==GOID,"name_1006" ][1]
		rm(GoTerm)
	}
	print(paste("Start Search", GOID, dateiName))	
	
	GoTermDF <- data.frame(getBM(attributes=c("entrezgene_id","ensembl_gene_id_version","external_gene_name","go_id"),
								 filters=c('go_parent_term'),
								 values = GOID ,
								 mart=ensembl))
	
	GoTermGeneList <- GoTermDF$external_gene_name 
	return(list(unique(GoTermGeneList),dateiName))
}

print(inputPath)
print(outputPath)

if (TRUE) {
	x <- scan(inputPath, what="", sep="\n")
	y <- strsplit(x, "[[:space:]]+")
	
	z <- unlist(y)
	listOfGO <- c()
	for (n in z) {
		if (substr(n,1,3)=="GO:") {
			listOfGO <- append(listOfGO, n)
		}
	}
	rm(x,y,z)
	
	goCollect <- c()
	for (n in listOfGO) {
		goGeneList <- GOTERMListFunction(n)
		goCollect <- append(goCollect, goGeneList[[1]] )
	}
	goCollect <- unique(goCollect)

	
	save(goCollect, file=outputPath)
}
