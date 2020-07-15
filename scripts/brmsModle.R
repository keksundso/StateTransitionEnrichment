#run the binomial modle
###############INPUTS and Outputs#########################
x = load(snakemake@input[["compDFName"]])
#x <- load("/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/ZS/threeGOList_completeDF.rda")
completeDF <- get(x) 
rm(x)

###############LIBRARIES#########################
library(brms)

#################CODE############################


iterValue <- 50000
warmupValue <- 5000

brmsModle <-
	brm(data = completeDF, family = binomial,
		hits | trials(trails) ~  (1 + GO | mut) ,
		iter = iterValue, warmup = warmupValue
		,  chains = 2, seed = 11, cores= 3
		, control = list(adapt_delta = 0.95,max_treedepth = 12)
	)

save(brmsModle, file= snakemake@output[["brmsModleVar"]])
write.csv2(brmsModle$prior, file = snakemake@output[["priorTable"]])