---
title: Codon models
subtitle: Substitution models for protein-coding sequences
authors:  Nicolas Lartillot
level: 3
order: 1
prerequisites:
- TP6_CodonModels
index: true
title-old: RB_CTMC_Tutorial
redirect: false
---

{% section Introduction %}

Thus far, the substitution models that we have considered are specified at the nucleotide level. They assume each nucleotide site evolves independently. In addition, they do not make the difference between synonymous and non-synonymous substitutions. By running these models on all coding positions or only on third coding positions, we have obtained indirect evidence suggesting that the substitution process is quite different between synonymous and non-synonymous changes. This is also expected based on our prior knowledge about evolutionary processes: synonymous substitutions are mostly reflecting the mutation rates, whereas non-synonymous substitutions incorporate the effect of selection at the amino-acid level, which is usually strong (thus, implying lower substitution rates) and also quite variable between sites.

In this tutorial, we will build a substitution model that is expressed directly at the level of codons. Thus, now, each position is a nucleotide triplet. There are 4x4x4 = 64 different triplets. However, stop codons are in principle excluded by very strong purifying selection in normal coding sequences, so we can exclude them. Under the universal genetic code, there are 3 stop codons. As a result, there are 61 non-stop codons. Altogether, the rate matrix $Q$ is now a 61x61 matrix, specifying the rate of substitution between each possible pair of non-stop codons.

There are several ways the codon substitution model could be specified. Here, we will use a simple approach, based on the following idea:

First, you specify a nucleotide mutation matrix $R$, which is 4x4. In the following we will use a HKY model for this. Thus, $R$ is equal to:

$$R = \begin{pmatrix}
* & \pi_{C} & \kappa \pi_{G} & \pi_{T} \\
\pi_{A} & *  & \pi_{G} & \kappa \pi_{T} \\
\kappa \pi_{A} & \pi_{C} & *  & \pi_{T} \\
\pi_{A} & \kappa \pi_{C} & \pi_{G} & *
\end{pmatrix} \mbox{  ,}$$

Now, consider two codons that differ only at one position. For instance, CCA and CCC. These two codons are synonymous: both code for Proline. Thus, the change from CCA to CCC corresponds to a single nucleotide change, from A to C, which is synonymous. Assuming it is neutral (no selection on codon usage, for instance), then it should occur at a rate equal to the mutation rate, which is here equal to $R_{AC}$. Here, $R$ is a HKY process, and A->C is a transversion, occurring at rate $\pi_{C}$ according to the matrix given above, so the rate from CCA to CCC is:

$$
Q_{CCA, CCC} = R_{AC} = \pi_C
$$

Similarly, the transition from CCA to CCG is also synonymous, but now it is a transition, which occurs, at the mutational level, at a different rate: $R_{AG} = $\kappa \pi_G$. Thus:

$$
Q_{CCA, CCG} = R_{AG} = \kappa \pi_G
$$

Now, consider the two codons CCA and CTA. The first encodes a proline, the second a leucine. The change is expected to be proposed by mutation at a rate $R_{CT} = \kappa \pi_T$. However, it is non-synonymous, and is thus probably under selection, so the rate of substitution (of change at the level of the whole population) is expected to differ from the mutation rate. If selection is purifying (if non-synonymous changes tend to be deleterious), then we expect that the substitution rate will be lower than the mutation rate. If, on the other hand, there is positive selection (e.g. in the context of a race between host and pathogens, new amino-acid variants might often be selected, for instance, if they lead to a more efficient escape from the defense of the host, etc), then the non-synonymous changes will more easily reach fixation in the population, and as a result, the substitution rate will be higher than the mutation rate.

To mathematically represent this, we invoke a new parameter, called $\omega$, which acts multiplicatively on non-synonymous substitutions. Thus, in the present case:

$$
Q_{CCA, CTA} = R_{CT} \omega = \kappa \pi_T \omega
$$

And we invoke this $\omega$ multiplicator for all non-synonymous substitutions. Thus:

$$
Q_{CCA, GCA} = R_{CG} \omega = \pi_G \omega
$$

And we do this for all possible single nucleotide changes.

Finally, we consider that double mutations are very unlikely. As a result, the rate of substitution, say, from AAA to ACC, or from AAA to CCC, is equal to 0.

Altogether, in order to build the 61x61 rate matrix $Q$, we have to consider all possible pairs of codons (C,D). For each of them, we determine whether they differ by only one nucleotide. If this is not the case, then we set $Q_{CD} = 0$. Otherwise, we determine the single nucleotide change that would lead from C to D (say, from nucleotide X to nucleotide Y), we read out the mutation rate for this nucleotide change from the 4x4 rate matrix $R$. Next, we determine whether the change is synonymous or non-synonymous. If synonymous, we say that $Q_{CD} = R_{XY}$. If non-synonymous, then $Q_{CD} = R_{XY} \omega$.

{% section Implementation %}

In order to run a codon analysis in RevBayes, we first need to read the nucleotide sequence alignment:
```
data <- readDiscreteCharacterData("data/placZFX20.nex")
```
Next, we should transform it into a codon sequence alignment. That is, we should explicitly say that each position is now a triplet of nucleotides:
```
data_codon = data.translateCharacters("Codon")
```
You can check that this indeed changes the way RevBayes considers the alignment by looking at the number of positions for the two alignments:
```
print(data)
print(data_codon)
```

In RevBayes, we could carefully construct the codon matrix, as described in the previous section. However, there is already a built-in function for this. 
```
Q := fnCodonHKY(omega=omega, kappa=kappa, baseFrequencies=nucstat)
```
This function, `fnCodonHKY`, directly takes the parameters of the nucleotide substitution model ($\kappa$ and $\pi$, in the mathematical notation above), as well as the $\omega$ parameter, and then automatically builds the 61x61 rate matrix $Q$.

Of course, before calling this function, you need to create the model variables `kappa`, `nucstat` and `omega`. For `kappa` and `nucstat`, you can do as previously, for the HKY model (an exponential prior for `kappa` and a uniform Dirichlet prior for `nucstat`). For `omega`, we can also use an exponential prior, which we could take of mean 5 (the dN/dS at the level of a gene is rarely above 1, in fact).

Once these three parameter components are defined, you can then create the matrix $Q$ (with the command given above), and the substitution process can then be instantiated:
```
seq ~ dnPhyloCTMC( tree=psi, Q=Q, type="Codon" )
seq.clamp( data_codon )
```
Note that we have to specify that the data are codons, not nucleotides (type = "Codon"). In addition, we have to clamp the substitution process to the codon-transformed data (`data_codon`) and not `data`, which corresponds to the original nucleotide alignment).

For the rest, the script is as before, except that we need to define a move on omega (a scaling move). 

Write the script, and run it on the ZFX gene (`placZFX20.nex`). Give a point estimate and a credible interval for omega. Do the same thing for BRCA1 (`placBRCA120.nex`). Again, estimate omega. How would you interpret the difference between ZFX and BRCA1?


{% subsection Increasing computational speed (optional) %}

As you can probably notice, the MCMC for this model is much slower than for the nucleotide-level analyses that we have conducted previously. The algorithmic complexity of the likelihood computation (i.e. the time it takes to compute the likelihood) is proportional to the square of the number of states: here, 61x61, as opposed to 4x4. There are three times less sites, and thus, in the end, the relative computational complexity between codon and nucleotide models is 61x61 versus 4x4x3, i.e. 3721 versus 192, which gives a factor 77.5. Each cycle takes 77.5 more time than under a nucleotide-level analysis!

There are other MCMC approaches that circumvent this problem, but they are not easily implemented or used in RevBayes. Here, one way to increase the speed of the inference is to fix the tree topology. Suppose that you have computed the MAP tree (or the consensus tree) for the ZFX dataset using a nucleotide model (e.g. T92). This tree is stored in a file, let us say that this file is called `ZFXmap.tree`. Then, what you can do is load the tree from file:
```
tree <- readTrees("analyses/ZFXmap.tree", treetype="non-clock")[1]
```
then, after creating the tree topology:
```
out_group = clade("Sorex")
topology ~ dnUniformTopology(taxa, outgroup = out_group)
```
you can constrain it (clamp it) to the topology specified by the MAP tree:
```
topology.clamp(tree)
```
Finally, you should deactivate the moves on the tree topology (i.e. remove the SPR and NNI moves). This will accelerate the MCMC, first, because less moves are done per cycle, but also, because the MCMC starts directly on a good tree.

Run the script thus modified on ZFX and on BRCA1.


{% section Modeling variation in omega across sites %}


In this section, we would like to model the fact that different coding sites have different selective constraints: some of them are highly constrained, while other are less so. As a result, we would expect to have different values of omega across sites.

A first approach to model this would be to do the same thing as for nucleotide models with rate variation among sites: use a discretized gamma distribution for omega. This model could be implemented in RevBayes. This is left as an optional exercise.

An alternative approach is to implement a mixture model. We specify three categories: one for purifying selection (with $\omega < 1$), one for neutral evolution ($\omega = 1$) and one with positive selection (with $\omega > 1$). As a consequence, we have three rate matrices, $Q_1$, $Q_2$ and $Q_3$, which are all parameterized by the same nucleotide-level parameters, but which differ in their $\omega$ parameter. Then, we specify mixture weights, by invoking a vector $w = (w_1, w_2, w_3)$, such that $w_1 + w_2 + w_3 = 1$. This vector belongs to the simplex S3. Finally, we say that sites across the protein are from this mixture:
- a proportion $w_1$ of them are under purifying selection ($dN/dS = \omega_1 < 1$);
- a proportion $w_2$ of them are under neutral evolution ($dN/dS = \omega_2 = 1$);
- a proportion $w_3$ of them are under positive selection($dN/dS = \omega_3 > 1$);

Once the model is specified, we can estimate the parameters, and in particular the proportions of sites under each of the three regimes, as well as the strength of purifying and positive selection ($\omega_1$ and $\omega_3$).

Ideally, we should also be able to estimate the posterior probability that each site belongs to each category. Of particular interest would be the posterior probability $p_i$ that site $I$ is under positive selection. However, the current implementation of RevBayes does not allow us to get this information.

{% subsection Implementation %}

We should first define our priors. For $\omega_1$, one could use a uniform prior between 0 and 1. The second entry, $\omega_2$, is fixed anyway. For $\omega_3$, one could define it as 1.0 + \delta$, where $\delta$ is a positive real number, which could have an exponential prior.

Then, we should create a vector of omega values, which we could call `omega_vector`, such that each of its entries would correspond to $\omega_i$ for $I=1, 2, 3$. However, we have to be careful: $\omega_1$ is a random variable, $\omega_2$ is a fixed constant, and $\omega_3$ is the sum between a fixed number and a random variable. In RevBayes, a vector should be entirely made of either random variables, or deterministic variables, so we can't define the entries of `omega_vector` directly. Instead, we need to first define the random variables separately, and then defined `omega_vector` as a vector of deterministic model variables.

Here, we first define two random variables, called `pur_om` (purifying selection, corresponding to $\omega_1$), and `delta_pos_om` (the variables called $\delta$ above, in the expression for $\omega_3$). Using the priors suggested above:
```
pur_om ~ dnBeta(1.0, 1.0)
delta_pos_om ~ dnExponential(1.0)
```

Next, we define `omega_vector` as a vector of deterministic variables:
```
omega_vector[1] := pur_om
omega_vector[2] := 1.0
omega_vector[3] := 1.0 + delta_pos_om
```

The weights of the mixture can be arbitrary. We can use a uniform Dirichlet prior (as we did for the nucleotide frequencies, but now with 3 entries, not 4):
```
omega_center = [1.0/3, 1.0/3, 1.0/3]
omega_weight ~ dnDirichlet(omega_center)
```

Next, one should create the 3 rate matrices, one for each value of omega contained in omega_vector:
```
for (i in 1:3)	{
	Q_vector[I] :=  fnCodonHKY(omega=omega_vector[i], kappa=kappa, baseFrequencies=nucstat)
}
```

Finally, we can specify the substitution process:
```
seq ~ dnPhyloCTMC( tree=psi, Q=Q_vector, siteMatrices = omega_weight, type="Codon" )
```
Several things are important here. First, the `Q` keyword is now associated, not with a single matrix, but with a vector of 3 matrices. Second, we have this additional argument `siteMatrices = omega_weight`, which specifies that the 3 matrices that were previously associated with the `Q` keyword should be understood as a mixture across sites, with weights given by `omega_weight`. When we do this, RevBayes will automatically average the likelihood at each site over the three components of the mixtures.

Once this is done, the model can be conditioned on the codon data, and MCMC moves should be defined for all of the free variables of the model.

Write the script and run it on ZFX and on BRCA1.

One question of particular interest is to ask whether there are sites under positive selection in a given gene. One way to address this question is to examine the posterior estimate of `omega_weight[3]`: this weight is our estimate of the proportion of sites for which $dN/dS > 1$. 

Would you say that ZFX has sites under positive selection? What about BRCA1 ? Would you say that purifying selection is stronger/weaker on the remaining sites, for ZFX, compared to BRCA1 ? 




