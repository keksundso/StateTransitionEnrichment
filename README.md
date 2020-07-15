
# StateTransitionEnrichment

To find state transitions which are enriched or depleted in given set of genes (compared to the Background), we analyzed our data via binomial regression in a bayesian framework.

p(T, GO) = aT + eT*GOT.

This way we are able to compare the probability that a given transition (pT) occurs in a gene which is part of a selected group (e.g. GO-Term) to all annotated genes. This difference is expressed be the effect size (e) and can be interpreted as the log odds ratio. In consequence this means that a negative log odds ratio shows us a depletion of any given transition in our selected group of genes, while a positiv value indicates an enrichment.


## Explaining the Analysis

Performing the analysis consists out of three branches ( see dag.svg).

1. Testing the accuracy of the model.
2. Testing the plausibility of the results.
3. Actual analysis of the genes of interested.

The accuracy of the model is assessed by posterior predictive checks (PPC) and the ability of the model to recover known, because simulated, data. The plausibly of the results is assessed by analyzing a set of genes (house keeping genes) where a set of transitions are expected to be enriched or depleted.

Explanation of the individual plots can be seen in the ***\*.rst*** files in **report/** folder.

## Prerequisite
### Software
Tool is build to be run in a  [conda enviroment](https://docs.conda.io/en/latest/) with a [snakemake pipline](https://snakemake.readthedocs.io/en/stable/index.html). The needed packages are defined in "condaPackages.txt". 

### Input-Data
* **ListOfGoTerms** a plain text file, format can be seen in example file **GOthreeComb.go**

* **journal.pone.0000898.s001.csv**, file in puplication  xxx can be retrived from and should be saved in Data folder is lifted to mouse genome, file can be requested. 

* **Jan2020gcsc3.mutations.RData**, file can be requested. 

## Executing the code with
    $ snakemake

