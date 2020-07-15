To analyze how well the model can predict the enrichment or depletion of a given transition, 334 genes of interested and 20000 genes of the background where simulated with a known probability of occurrence for 15 different states (a-o) in two cell states.
The binomial regression model was trained with this dataset and calculated the posterior distribution of log odds ratios (effect sizes) for each of the 225 transitions.

Here the most likely effect sizes are shown on the ordinate, as well as the 11% — 89% density intervals of the posterior distribution.
On the abscissa the difference of the simulated probabilities logit values are shown, which correspond to the log odds ratio for the given transition.
It can be seen that there is a high correlation between the simulated and the recovered data, inferring a high predictive accuracy of the model. Only very large differences (larger than the ones observed in the actual analysis) show a meaningful divergence from the diagonal.