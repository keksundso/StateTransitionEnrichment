
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('ZS/threeGOList_completeDF.rda', "compDFName" = 'ZS/threeGOList_completeDF.rda'),
    output = list('ZS/threeGOList_BinomialModle.rda', "brmsModleVar" = 'ZS/threeGOList_BinomialModle.rda'),
    params = list(),
    wildcards = list('threeGOList', "GeneListWK" = 'threeGOList'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list(),
    rule = 'BRMSModle'
)
######## Original script #########
x = load(snakemake@input[["compDFName"]])
#x <- load("/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/ZS/threeGOList_completeDF.rda")
raw <- get(x) 
rm(x)

library(brms)

iterValue <- 1000
warmupValue <- 200

brmsModle <-
	brm(data = completeDF, family = binomial,
		hits | trials(trails) ~  (1 + GO | mut) ,
		iter = iterValue, warmup = warmupValue,  chains = 2, seed = 10, cores= 3
		, control = list(adapt_delta = 0.95,max_treedepth = 12)
	)

save(brmsModle, file= snakemake@output[["brmsModleVar"]])