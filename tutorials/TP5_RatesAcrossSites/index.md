---
title: TP5: Accounting for variation in substitution rates across sites
subtitle: A mixture model to account for rate variation across sites
authors:  Nicolas Lartillot
level: 2
order: 0.5
prerequisites:
- TP5_RatesAcrossSites
index: true
title-old: RB_CTMC_Tutorial
redirect: false
---


Thus far, all phylogenetic models that we have considered assume that all sites evolve at the same rate and under the same substitution process (such as specified by the rate matrix $Q$). In reality, different sites are under different constraints. As a result, some sites might be fast evolving, while other sites might be slowly evolving. We need to account for this variation across sites in our models.

{% section A random-effect model with gamma-distributed rates across sites %}

The most widely used approach to model varying rates across sites is to use a random-effect model. The idea is to assume that each site has a certain rate of evolution (let us call it $r_i$, for site $i$). These rates are assumed to be independent and identically distributed across sites from a parametric distribution. Generally, we assume a gamma distribution.

A gamma distribution has two parameters: a shape parameter ($\alpha$) and a scale parameter ($\beta$). You can have a look at the Wikipedia page about the gamma distribution to have more information. The mean of the distribution is equal to $\alpha / \beta$, and the variance is $\alpha / \beta^2$. Thus, by tuning the values of these two parameters, $\alpha$ and $\beta$, we can create a broad variety of distributions, with arbitrary mean and variance.

Here, we want to model the relative rates across sites (we don't need to model the absolute rates because we have the branch lengths, or the total tree length, which will capture the overall rate of evolution). So, we can constrain the mean to be equal to 1, which is equivalent to setting $\beta = \alpha$. As a result, our distribution has only one parameter, $\alpha$, and the variance of the distribution is equal to $\alpha / \alpha^2 = 1 /\alpha$.

Computationally speaking, averaging the likelihood over a gamma distribution at each site is difficult. For that reason, most often, the gamma distribution is discretized, into a small number of categories (typically, as few as $K=4$ categories). As a result, the integral at each site can be done by simply averaging the likelihood over the $K$ categories.


{% section Implementing the model in RevBayes %}

We don't know the value of the shape parameter $\alpha$, or equivalently, we don't know the variance of the gamma distribution (which, as we said, is the inverse of $\alpha$). So, we need to consider it as a parameter, with a prior. Here, we parameterize directly in terms of the variance. We use a broad exponential prior, of mean 5 (thus, of rate $\lambda = 0.2$):
```
rate_var ~ dnExponential(0.2)
```
Then, we set $\alpha$ as a deterministic variable, equal to the inverse of the variance:
```
alpha := 1.0 / rate_var
```

Next, we create a discretized gamma distribution with 4 categories, of shape and scale parameters both equal to $\alpha$:
```
disc_rates := fnDiscretizeGamma(alpha, alpha, 4)
```

Finally, we specify the substitution process. Now, in addition to the tree and the substitution rate matrix $Q$, we also specify the distribution rates across sites with the 'siteRates' keyword:
```
seq ~ dnPhyloCTMC(tree=psi, Q=Q, siteRates=disc_rates, type="DNA" )
```

For the rest, the script is as before. The rate matrix Q could be anything: either the JC, or the T92 or the HKY or the GTR model. Here, we will analyse the BRCA1 gene (again, and short version), for which we can use the T92 model. Concerning the moves, don't forget to add a move (a scaling move) for the 'rate_var' parameter.

Draw the graph of this model. Write the script and run it on the BRCA1 gene. Try first with all coding positions ('BRCA1short.nex'), and then with only the third coding positions ('BRCA1third.nex'). Run for 10 000 generations. Using Tracer, visualize the traces, determine the burnin, visualize the posterior distribution on the variance parameter and obtain a point estimate and a 95\% credible interval. Do this for the two datasets (all coding positions and third coding positions). How do you interpret the difference in the variance parameter between these two datasets?

Do you see major differences in tree inference, compared to the model without rates across sites?


