******************************
--- Explaining the model
******************************

To find state transitions which are enriched or depleted in given set of genes (compared to the Background), we analyzed our data via binomial regression in a bayesian framework.

:math:`p_{(T,GO)} = a_T + e_T * GO_T`.

This way we are able to compare the probability that a given transition (:math:`p_{T}`) occurs in a gene which is part of a selected group (e.g. GO-Term) to all annotated genes.
This difference is expressed be the effect size (:math:`e`) and can be interpreted as the log odds ratio.
In consequence this means that a negative log odds ratio shows us a depletion of any given transition in our selected group of genes, while a positiv value indicates an enrichment.

******************************
--- Explaining the Analysis
******************************

Performing the analysis consists out of three branches (dag.svg_ ).

1. Testing the accuracy of the model.
2. Testing the plausibility of the results.
3. Actual analysis of the genes of interested.

The accuracy of the model is assessed by posterior predictive checks (PPC) and the ability of the model to recover known, because simulated, data.
The plausibly of the results is assessed by analyzing a set of genes (house keeping genes) where a set of transitions are expected to be enriched or depleted.