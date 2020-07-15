# creating a list of houskeeping genes
#https://doi.org/10.1371/journal.pone.0000898 this is the source for the house keeping genes



###############INPUTS and Outputs#########################
HKGTextPath <- snakemake@input[["HKGcsv"]]
#HKGTextPath <- "/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/data/ListOfGenes/journal.pone.0000898.s001.csv"

###############LIBRARIES#########################
#integrated in the code section 

#################CODE############################
convertHumanGeneList <- function(x){
	
	require("biomaRt")
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	
	genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	
	humanx <- genesV2[, 2]
	
	# Print the first 6 genes found to the screen
	print(head(humanx))
	return(genesV2)
}



#https://doi.org/10.1371/journal.pone.0000898
library(readr)
HK_df <- data.frame(read_delim("/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/FremdDaten/HKG/journal.pone.0000898.s001.csv",  "\t", escape_double = FALSE, trim_ws = TRUE))
HK_df <- HK_df[HK_df$MFC < 2 & HK_df$CV < 10, ]

rownames(HK_df) <- HK_df$Identifiers

mouse_genes <- convertHumanGeneList(HK_df$Identifiers)

HK_df <- merge(mouse_genes, HK_df, by.x = 1, by.y = 1)


##############BioMart annotation#######################
library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#listDatasets(ensembl)
attributes = listAttributes(ensembl)
#attributes[1:100,]
filters = listFilters(ensembl)
#filters[1:100,]

N <- data.frame(getBM(attributes=c("mgi_symbol","ensembl_gene_id_version","external_gene_name","chromosome_name","strand"),
					  filters=c('mgi_symbol'),
					  values = HK_df$MGI.symbol,
					  mart=ensembl))

names(HK_df)[names(HK_df) == 'MGI.symbol'] <- 'mgi_symbol'
HK_df <- merge(N,HK_df,by="mgi_symbol")
#save(K_df, file = "KallistoOutputAnnotatet.RData")
rm(attributes,ensembl,filters,N,mouse_genes)


ListOfHKG <- HK_df$external_gene_name

save(ListOfHKG, file=snakemake@output[["HGKVector"]])