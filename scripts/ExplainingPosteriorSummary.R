# create an explanation plot, which describes how to each data point in the transition enrichment plot came to be
###############INPUTS and Outputs#########################

testRun <- FALSE

if (testRun == TRUE) {
 #not needed jet.
}
if (testRun == FALSE) {
	x = snakemake@input[["brmsModleVar"]]
	
	outPutPath = snakemake@output[["plotExp"]]
}

brmsModle <- get(load(x))
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

pdf(outPutPath)

raw <- posterior_samples(brmsModle)$`r_mut[10.9,GO]`

x.norm <- seq(-2, 2, by = .001)
y.norm <- dnorm(x.norm, mean = mean(raw), sd = sd(raw))


plot(density(raw), lty = 3, col= "gray", main = "14.10",xlab= "log(Odds-Ratio)")
lines(x.norm,y.norm)

lines(c(myMod(raw),myMod(raw)), c(0.1,10))
text(myMod(raw), 0.05, labels = round(myMod(raw), digits = 2) )

polygon(c( x.norm[x.norm<=0], 0 ),  c(y.norm[x.norm<=0],0 ), col="gray")
text(0.0, 0.05, labels = round(pnorm(0, mean =mean(raw), sd = sd(raw), lower.tail = TRUE, log.p = FALSE), digits = 2), pos=2)

legend("topleft", legend = c("Empirical Distribution", "Normal Fit"), lty = c(3,1), col= c("gray","black"))

dev.off()
