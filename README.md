
# StateTransitionEnrichment

To find state transitions which are enriched or depleted in given set of genes (compared to the background), we analyzed our data via binomial regression in a bayesian framework.

Yt,x ~ binomial(θt,x , Nx) 
Θt = logistic(αt+βt*Xt) 

The number of genes with transition t and the condition x (part of the background (x = 0) or the selected group (x = 1)) is Yt,x. It depends on the number of genes in condition x (Nx) and the probability θt,x that a given transition t occurs in a gene which is annotated to one of the two conditions x. The influence of the condition on θ is expressed by the effect size β, which is the log odds ratio. The uncertainty that an enrichment (positive log odds ratio) or a depletion (negative log odds ratio) occurs is represented by the probability that β is either negative or positive, respectively.


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

