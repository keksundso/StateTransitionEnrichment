#Show how well the simulated data gets recovered by the model.
###############INPUTS and Outputs#########################

testRun <- FALSE

if (testRun == TRUE) {

	simModle <- "/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/ZS/Sim_BinomialModle.rda"
	
	propDF <- "/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/ZS/SimPropDF.rda"
	
	
	outPutPathPDF <- "/Volumes/Entwicklungsbiologie/Christoph/PhD-Work_DataUScrips/H_vs_P_RNASeq_ChromHMM/h_vs_p_rnaseq_chromhmm_2019/mutationGeneMap/baysianApproch/Workflow_snakemake/output/SimRecovery.pdf"
}
if (testRun == FALSE) {
	simModle = snakemake@input[["SimModlePath"]]
	propDF = snakemake@input[["propDFVar"]]
	outPutPathPDF = snakemake@output[["simRecPath"]]
}

propDF  <- get(load(propDF))
meinModel <- get(load(simModle))

###############LIBRARIES#########################

library(brms)

#################CODE############################





plotPredRecov <- function(dasModel){
	# Create PPC data table ready to be plotted
	
	postSampleDF <- posterior_samples(dasModel)
	listofx <- c()
	listofy <- c()
	counter  <- 0
	y1 <- c()
	y2 <- c()
	flag <- c()
	for (mut in ListOfMutations) {
		conter <- counter +1
		
		for(n in colnames(postSampleDF)[grepl(mut,colnames(postSampleDF)) & grepl("GO",colnames(postSampleDF))&grepl("r_mut\\[",colnames(postSampleDF))]){
			counter <- counter + 1
			y <-mean(postSampleDF[,n])
			
			if (TRUE) { # 11 and 89% quantile  version
				tmp <-  quantile((postSampleDF[,n]), probs = c(0.11,0.89))
				#print(paste(mut, tmp[1],mean(postSampleDF[,n]),tmp[2]), sep= "   --   ")
				y1 <- c(y1,   tmp[1]  )
				y2 <- c(y2,   tmp[2]  )
			}
			
			if (FALSE) { # sd version
				y1 <- c(y1,   mean(postSampleDF[,n])- sd(postSampleDF[,n])   )
				y2 <- c(y2,   mean(postSampleDF[,n])+ sd(postSampleDF[,n])   )
			}
			
			
			
			
			x <- logit_scaled(propDF[mut, "GOProp"]) - logit_scaled(propDF[mut, "BGProp"] )
			
			if (propDF[mut, "GOProp"] == 0 | propDF[mut, "GOProp"] == 1) {flag <- c(flag, "red")	}
			if (propDF[mut, "GOProp"] != 0 & propDF[mut, "GOProp"] != 1) {flag <- c(flag, "black")	}
			
			#plot(density(postSampleDF[,n]))
			listofx <- c(listofx, x)
			listofy <- c(listofy, y)
			
			#if (x < -4 | x > 4) {print(paste(mut, n))	}
		}
	}
	
	
	return(list(listofx,listofy,y1,y2))
}


ListOfMutations <- c()
for (n in letters[1:15]) {
	for (m in letters[1:15]) {
		ListOfMutations <- c(ListOfMutations, paste0(n,".",m))
	}
}

if (testRun == FALSE) {pdf(outPutPathPDF,width=8,height=8)}

Uberschrift <- basename(simModle)
if (TRUE) {# create the plot
	functionOut <- plotPredRecov(meinModel)
	par(pty="s")
	
	maxValue <- max(abs(unlist(functionOut)))
	
	plot(functionOut[[1]],functionOut[[2]], ylab= "log(ODDs-Ratio) (Recovered)", xlab="log(ODDs-Ratio) (Simulated)", pch="."
		 , xlim= c(-maxValue,maxValue)
		 , ylim= c(-maxValue,maxValue)
		 ,main=Uberschrift
	)
	counter <- 0
	for (n in functionOut[[1]]) {
		counter <- counter +1
		#print(c(x1[counter],x2[counter]))
		lines(c(functionOut[[1]][counter],functionOut[[1]][counter]),c(functionOut[[3]][counter],functionOut[[4]][counter]))
		
	}
	text(functionOut[[1]],functionOut[[2]],ListOfMutations, cex=0.5)
	lines(c(100,-100),c(100,-100), lty=3)
}


if (testRun == FALSE) {dev.off()}