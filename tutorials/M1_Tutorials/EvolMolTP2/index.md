---
title: Evol Mol TP2 - Nucleotide substitution models
subtitle: Estimating substitution rates between nucleotides
authors:  Nicolas Lartillot
level: 3
order: 0.4
prerequisites:
- EvolMolTP1
index: true
title-old: RB_CTMC_Tutorial
redirect: false
---

The last tutorial was introducing how to conduct a Bayesian phylogenetic analysis under the Jukes-Cantor model. In this tutorial, we will follow up on this, by exploring other models of nucleotide substitution.

A model of nucleotide substitution specifies the rates of substitution between all pairs of nucleotides. It is entirely defined by a 4x4 matrix called instantaneous rate matrix, or more simply, rate matrix ($Q$).

We have seen in the last tutorial that the instantaneous rate matrix for the JC model is given by:

$$Q_{JC} = \begin{pmatrix}
{*} & \frac{1}{3} & \frac{1}{3} & \frac{1}{3} \\
\frac{1}{3} & {*} & \frac{1}{3} & \frac{1}{3} \\
\frac{1}{3} & \frac{1}{3} & {*} & \frac{1}{3} \\
\frac{1}{3} & \frac{1}{3} & \frac{1}{3} & {*}
\end{pmatrix} \mbox{  ,}$$

This model, however, is simplistic. In practice, different types of substitutions (or different types of point mutations) might not occur at the same rate. Thus, we might need to use more sophisticated models.

Here, we will explore 3 alternative models, of increasing complexity. For each substitution model, we will write a script and run it. This will give us a joint posterior distribution on trees, but also on model parameters. Once the analysis has been conducted, we can then interpret the posterior distribution over the parameters of the model. This will give us information about the process of molecular evolution.

The models explored in this tutorial are summarized in this table.

{% table tab_subst_models %}
{% tabcaption %}
Specific functions for substitution models available in RevBayes.
{% endtabcaption %}

 |   **Model**      |        **Reference**        |  **Function**   |      **Parameters**     |
 |:----------------:|:---------------------------:|:---------------:|:-----------------------:|
 |   Jukes-Cantor   |    {% cite Jukes1969 %}     |      fnJC       |           -             |
 |        T92       |    {% cite Tamura1992 %}    |      fnT92      | $\pi_{GC}$, $\kappa$    |
 |        HKY       |   {% cite Hasegawa1985 %}   |      fnHKY      |     $\pi$, $\kappa$     |
 |        GTR       |    {% cite Tavare1986 %}    |      fnGTR      |     $\pi$, $\epsilon$   |

{% endtable %}
They are now presented in more details.

{% section The Tamura 92 model (T92) %}

{% subsection Description of the model %}

This model has two parameters, $\kappa$ and $\gamma$. The rate matrix is given by:

$$Q_{T92} = \begin{pmatrix}
{*} & \frac{\gamma}{2} & \kappa \frac{\gamma}{2} & \frac{1-\gamma}{2} \\
\frac{1-\gamma}{2} & {*} & \frac{\gamma}{2} & \kappa \frac{1-\gamma}{2} \\
\kappa \frac{1-\gamma}{2} & \frac{\gamma}{2} & {*} & \frac{1-\gamma}{2} \\
\frac{1-\gamma}{2} & \kappa \frac{\gamma}{2} & \frac{\gamma}{2} & {*}
\end{pmatrix} \mbox{  ,}$$

Two things are important to note:
- if $\kappa > 1$, then this says that transitions between purines (A and G) or between pyrimidines (C and T) are more frequent than transversions (from a purine to a pyrimidine, or vice versa). If $kappa < 1$, on the other hand, then transitions are less frequent than transversions.
- if $\gamma < 0.5$, then substitutions from GC to AT occur at a higher rate, compared to substitutions from AT to GC, and conversely if $\gamma > 0.5$.

For which values of $\kappa$ and $\gamma$ does the T92 model reduce to JC ?

{% subsection Implementation %}

A script (`scripts/phyloJC.rev`), implementing phylogenetic inference with the JC model, is given with the tutorial. It is slightly different from the one that you used for the last tutorial, but is essentially equivalent. It was just re-written so as to better separate the model declaration from the moves and the monitors. In addition, the move schedule was slightly modified, so as to allow for faster analysis on the specific dataset of interest for today. This dataset is an alignment of the BRCA1 gene in placental mammals (`data/placBRCA1.nex`). To further accelerate the analyses, we will in fact use a short version of this dataset (in the file called `placBRCA1short.nex`), which contains only the 600 first nucleotide positions of this gene. Finally, a third dataset is provided (file called `placBRCA1third.nex`), which contains only the third coding positions of the original alignment. In the following, we will explore both `placBRCA1short.nex` and `placBRCA1third.nex`.


Copy the JC script into a new script, which you may call `phyloT92.rev`), and open it. You will modify this script so as to implement a T92 model. Much of the script should stay the same, but you will have to introduce several modifications in order to take into account the specific aspects of this new and more complex model.

First, under JC, the rate matrix and the sequence evolutionary process were specified as follows:
```
Q <- fnJC(4)
```
In the case of JC, the rate matrix is a constant, hence the use of `<-`.

Then, we create a stochastic variable, which represents the multiple sequence alignment, such as produced by a substitution process running along the tree psi, and with substitution rates given by $Q$:
```
seq ~ dnPhyloCTMC(tree=psi, Q=Q, type="DNA")
```

If we now want to use a T92 model instead of a JC model, we need, first, to defined the two parameters: $\kappa$ and $\gamma$, as stochastic variables having a prior. Then, we can build the $Q$ matrix, and finally, we create a `phyloCTMC` (`seq`) using this $Q$ matrix. Finally, we need to move these two parameters during the MCMC.

The transition-transversion rate ratio $\kappa$ is a non negative real number. A reasonable prior assumption about $\kappa$ is that it is not much higher than 10. So, we can use an exponential prior of mean 10 (of rate $\lambda = 0.1$).
```
kappa ~ dnExponential(lambda = 0.1)
```
Concerning $\gamma$, it is a number between 0 and 1. We can invoke a uniform prior distribution. The uniform distribution between 0 and 1 is a particular case of the beta distribution, with with parameters 1.0 and 1.0:
```
gamma ~ dnBeta(1.0, 1.0)
```
Then, we can create the rate matrix $Q$, which is now a T92 rate matrix:
```
Q := fnT92(kappa = kappa, gc = gamma)
```
Note that the rate matrix is now a deterministic model variable (and not anymore a constant) -- hence, we use `:=`.

Finally, we can proceed with the sequence evolutionary process:
```
seq ~ dnPhyloCTMC(tree=psi, Q=Q, type="DNA")
```

Concerning the MCMC moves, we need to add moves for $\kappa$ and $\gamma$. For $\kappa$, we can use a scaling (or multiplicative) move, which is the best move to use for positive real numbers. For $\gamma$, we can use a sliding (or additive) move. As in the case of the JC model, we also need to move the tree topology, using NNI and SPR moves, as well as the branch lengths (using scaling moves). Altogether, or move vector can be specified as follows:
```
moves = VectorMoves()
moves.append(mvNNI(topology, weight=3.0))
moves.append(mvSPR(topology, weight=3.0))
moves.append(mvSlide(gamma, weight=1.0))
moves.append(mvScale(kappa, weight=1.0))
for (i in 1:n_branches) {
   moves.append(mvScale(bl[i], weight=1.0))
}
```

Write the script, run it on the BRCA1 dataset (short version). Load the log file in Tracer. Estimate the burnin and visualize the posterior distribution for $\kappa$. Give a point estimate (median) and a symmetric 95\% credible interval. Do the same thing for $\gamma$.

Based on your point estimated and credible intervals, would you say that T92 provides a better fit to the data, compared to JC?

Do a similar analysis, now with only the third coding positions (file called `placBRCA1third.nex`). Run the chain for 10 000 generations. Visualize the posterior distribution and give a point estimate and a credible interval for $\kappa$ and $\gamma$. Compare with the estimates obtained under the complete dataset (all coding positions). How do you interpret the differences?

Print out the map or posterior consensus tree. Does that differ a lot from the consensus tree inferred using the JC model?

{% section The Hasegawa, Kishino and Yano model (HKY)%}


The HKY model is slightly more complex than T92. It also relies on a transition/transversion ratio ($\kappa$). However, whereas T92 assumes that the equilibrium frequencies of G and C are the same (both equal to $\gamma/2$), and similarly for A and T (with frequencies both equal to $(1-\gamma)/2$, the HKY does not make this assumption. Under HKY, the equilibrium frequencies of the 4 nucleotides can be arbitrary. They are mathematically formalized by a frequency vector of dimension 4, $\pi = (\pi_A, \pi_C, \pi_G, \pi_T)$ such that:
$$\pi_A + \pi_C + \pi_G + \pi_T = 1$$

The rate matrix of HKY is:

$$Q_{HKY} = \begin{pmatrix}
{*} & \pi_{C} & \kappa \pi_{G} & \pi_{T} \\
\pi_{A} & {*}  & \pi_{G} & \kappa \pi_{T} \\
\kappa \pi_{A} & \pi_{C} & {*}  & \pi_{T} \\
\pi_{A} & \kappa \pi_{C} & \pi_{G} & {*}
\end{pmatrix} \mbox{  ,}$$


The T92 model is a particular case of the HKY model. For which values of its parameters does HKY reduce to T92?

The domain of definition of frequency vectors is called the simplex (here, the simplex S4, since there are 4 entries for the vector). A standard distribution on the simplex is the Dirichlet distribution (you can look at the definition of the Dirichlet distribution on Wikipedia). In particular, a uniform distribution over the simplex is a Dirichlet with weights all equal to 1. Thus, we can assume this prior for $\pi$:
```
pi ~ dnDirichlet([1.0, 1.0, 1.0, 1.0])
```

As for the T92 model, we can assume a broad exponential prior for $\kappa$:
```
kappa ~ dnexponential(lambda = 0.1)
```
and then define the HKY rate matrix:
```
Q := fnHKY(kappa=kappa, baseFrequencies=pi)
```
and finally, create the substitution process:
```
seq ~ dnPhyloCTMC(tree=psi, Q=Q, type="DNA")
```

We also need to move $\kappa$ and $\pi$. For $\kappa$, we can use a scaling move, as we did for the T92 model. For $\pi$, we can use a move that preserves the positivity of the entries of the vector, but also the constraint that the 4 entries of the vector should sum to 1. The move that does this is a `DirichletSimplex` move. The syntax for this move would be as follows:
```
moves.append(mvDirichletSimplex(pi, weight=1.0, alpha=10))
```

Write and run the script on the BRCA1 gene (only the third coding positions). Again, visualize the posterior distribution over the model parameters. In particular, visualize simultaneously the posterior distributions over the 4 nucleotide frequencies. Would you say that HKY should be used instead of T92? Is tree inference fundamentally different between the two models?

{% section The GTR model (optional) %}

The most general time reversible model can be parameterized as follows:
$$Q = \begin{pmatrix}
{*} & \rho_{AC} \pi_C & \rho_{AG} \pi_G & \rho_{AT} \pi_T \\
\rho_{CA} \pi_A & {*}  & \rho_{CG} \pi_G & \rho_{CT} \pi_T \\
\rho_{GA} \pi_A & \rho_{GC}  \pi_C & {*}  & \rho_{GT} \pi_T \\
\rho_{TA} \pi_A & \rho_{TC} \pi_C & \rho_{TG} \pi_G & {*}
\end{pmatrix} \mbox{  ,}$$
with the constraint that $\rho_{XY} = \rho_{YX}$ for all nucleotide pairs $(X,Y)$. The $\rho$'s are called the relative exchangeabilities (or exchange rates), and the $\pi$'s are the equilibrium frequencies of the process. With the constraint of symmetry for $\rho$, we have 6 distinct values (6 unordered pairs of nucleotides). In addition, and without loss of generality, we can constrain the $\rho$ vector to sum to 1. With these constraints, $\rho$ is a 6-dimensional frequency vector (i.e. it is a point on the simplex $S6$).

We can use a uniform prior over $\rho$, which is a particular case of the Dirichlet distribution with equal weights:
```
rho ~ dnDirichlet([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
```
and a uniform prior over $\pi$ (as for the HKY model):
```
pi ~ dnDirichlet([1.0, 1.0, 1.0, 1.0])
```
We can then define the GTR rate matrix:
```
Q := fnGTR(exchangeRates=rho, baseFrequencies=pi)
```
and create the substitution process:
```
seq ~ dnPhyloCTMC(tree=psi, Q=Q, type="DNA")
```
Finally, we need to move $\rho$ and $\pi$, in both cases using a `DirichletSimplex` move (as we did above for pi in the case of the HKY model).

Write the script and run it on the BRCA1 dataset (only the 3 coding positions). Explore the posterior distribution over the parameters and, in particular, determine whether all transitions occur at the same relative rate. Same thing for transitions.
