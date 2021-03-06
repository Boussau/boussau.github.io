---
title: TP2 Implementing Bayesian inference with RevBayes
subtitle: Inferring the divergence time between Human and Chimpanzee
authors:  Nicolas Lartillot
level: 2
order: 0.2
index: true
prerequisites:
- TP1_simulatingDNAEvolution
redirect: true
---


In this tutorial, we will conduct our first Bayesian analysis with RevBayes. Our aim will be to infer the divergence time between Humans and Chimpanzees.
To do this, we will come back to the simulation model that was developed in the previous tutorial.
This simulation model was taking the branch length and the substitution rate as the input, and was giving as an output the number of substitution events, or the number of differences between two taxa in their nucleotide sequences.

Here, we consider the problem in the other direction:
we have empirical data, in the form of a sequence alignment between two species (Human and Chimp). We can count the number of differences between the sequences of these two species. We also have some knowledge about the mutation rate in apes. Based on these data, we would like to infer the divergence time between Human and Chimp.



The data
===============
{:.section}


The file apes.nex contains an alignment with $N=878$ aligned nucleotide positions, between Human, Chimp and Gorilla. This alignment is a concatenation of 17 nuclear genes, for which only the 4-fold degenerate third codon positions have been kept. To visualize this alignment, you can open it with SeaView.
Out of these $N$ nucleotide positions, $k_{obs} =11$ are different between the human and the chimpanzee sequences.

We use only the 4-fold degenerate third codon positions because all mutations at these positions are synonymous, and therefore, we can assume that they are neutral (they are not under selection). As a result, the rate at which substitutions will accumulate over time is equal to the mutation rate.

Based on pedigree sequencing, the mutation rate per generation in Human has been estimated at $u = 2.10^{-8}$. Similar experiments have been conducted in Chimpanzees, giving similar estimates. The generation time in Humans is around $\tau = 30$ years. It is a bit shorter in Chimpanzees, however, for simplicity, we will ignore this, and we will assume that the generation time has always been 30 years in great apes.

Based on these estimates for $u$ and $\tau$, What is the mutation rate per million years? Let us call it $r$.


The Jukes-Cantor distance
===============
{:.section}


We would like to estimate the divergence time between Human and Chimp. Let us call this divergence time $T$ (measured in million years).

- Given $r$ and $T$, and based on what you did in the first tutorial, what is the probability that a given nucleotide position will be different between Human and Chimp?
Let us note this probability $q(r,T)$. It is a function of the rate $r$ and the divergence time $T$.

- Given a sequence of length $N$, what is the expected number $\bar k$ of differences between Human and Chimp over this sequence?

- Conversely, given that you have observed $k_{obs}$ out of $N$ differences in the empirical data, what is your best estimate for the divergence time $T$?

Of note, $\bar k$ is the expected number of differences, that is, the mean number of differences that you observe if you run the simulation many times under a given fixed value for $T$. However, from simulation to simulation, the actual number of differences $k$ will vary. This number $k$ is a random variable, which has a binomial distribution:
$$
k \sim Binomial(N, q(r,T))
$$
We can see this by noting that, for each position, we will end up observing a difference with probability $q(r,T)$. We repeat this procedure $N$ times independently, and therefore, the total number of 'successes' (i.e. differences) will be binomial of parameter $N$ and $q(r,T)$. This binomial distribution for $k$ is very important for the following.


Bayesian inference of the divergence time $T$: rejection sampling.
===============
{:.section}

The quick estimation obtained in the previous section does not give any measure of the uncertainty about our estimate of the divergence time $T$. To address this problem, we will now develop a more elaborate approach, based on Bayesian inference. This will also give us an occasion to write our first Bayesian analysis in RevBayes.

To conduct Bayesian inference, we first need to have a prior distribution over $T$. This prior represents our state of knowledge about $T$, before we have seen the empirical data (the sequence alignment). In the present case, we will consider that we don't know anything a priori about this divergence time, except that it is between 0 and 20 Myr. Mathematically, we represent this by assuming a uniform distribution over $T$. In RevBayes, we can draw $T$ from this prior as follows:
$$
T ~ Uniform(0,20)
$$
This is the mathematical notation. In RevBayes, we can draw $T$ from this prior as follows:

```
T = rUniform(1,0,20)[1]
```
You can try this command and repeat it several times, to see that, each time, you draw a value of $T$ anywhere between 0 and 20 Myr.

Now, consider the following program:

```
# number of aligned nucleotide positions
N = 878
# number of differences observed in the empirical sequences
k_obs = 11

# give your formula for r
r = ...

# number of samples we want from the posterior distribution
n_sample = 100

# initialize the counter
count = 0

# create a vector of size n_sample, with all entries initialized at 0
T_sample = rep(0, n_sample)

# repeat until we have collected n_sample values of T
while (count < n_sample)	{

	# draw T from prior: uniform between 0 and 20 Myr
	T = rUniform(1,0,20)[1]

	# compute q(r,T) for the current value of T
	# write some code here based on the formula that was derived previously
	q = ...

	# draw k given T
	k = rbinomial(1,N, Probability(q))[1]

	# condition on k == k_obs (thus rejecting the sample if k is not equal to k_obs)
	if (k == k_obs)	{
		count = count + 1
		T_sample[count] = T
	}
}

# finally, write the sample to a file:
for (i in 1:n_sample)	{
	write(T_sample[i], "\n", file="analyses/apes_reject.log", append=TRUE)
}
```

This program repeats the following procedure:
- draw T from the uniform prior
- given T, draw the number of differences k between the human and chimp sequences
- if k == k_obs, keep T, otherwise, discard it
- do this until 100 values for T have been kept

Import this program in RevBayes and run it. In R, draw the histogram of the 100 values of $T$ that have been stored in the array named sample (and that have been written into the file named apes_reject.out). To get a nicer histogram, you can set n_sample to 1000, but this will take more time.

In RevBayes, you can also compute the mean and the stdev of the sample:
```
mean(T_sample)
stdev(T_sample)
```
or the median
```
median(T_sample)
```
and, finally, have a 95\% credible interval by looking at the quantiles at 2,5\% and 97,5\%:
```
quantile(T_sample, 0.025)
quantile(T_sample, 0.975)
```

As you can see, although the values for $T$ that were initially drawn are uniformly distributed between 0 and 20, the histogram is concentrated around a value of about 9 Myr. In other words, this histogram is enriched in values of $T$ that are more likely to produce a value for the random variable $k$ which turns out to be equal to the empirically observed value $k_{obs}$.

From a Bayesian perspective, this histogram represents a sample from the posterior distribution. The algorithm used here to sample from the posterior distribution is called rejection sampling. Rejection sampling represents the concept of conditioning the model (including the prior) on the observed value $k_{obs}$ (simply by rejecting all runs that did not produce this result).


Bayesian inference of the divergence time $T$. MCMC.
===============
{:.section}


Rejection sampling is conceptually simple, but it is not a very efficient algorithm. For large datasets and more complex models, the probability of reproducing exactly the empirically observed data is just too small. For that reason, we will need to make a huge number of trials until we get 100 or 1000 samples.

A more efficient approach is to use Markov Chain Monte Carlo (MCMC). RevBayes semi-automatically implements MCMC. In this section, we will derive a MCMC version of the inference conducted above using rejection sampling.


Declaring the model
---------------------------------------------------
{:.subsubsection}


As above, we start by defining the fixed parameters (the constant) of the problem: $N$, $r$ and $k_{obs}$.

```
# number of aligned nucleotide positions
N <- 878
# number of differences observed in the empirical sequences
k_obs <- 11

# give your formula for r
r = ...
```


Then, we specify the model:

```
# T is from a uniform prior
T ~ dnUniform(0,20)

# here, your formula for q(r,T)
q := Probability(...)

# k is binomial, given T (q depends on T)
k ~ dnBinomial(N, q)

# k has been observed and turns out to be equal to k_obs
k.clamp(k_obs)

mymodel = model(T)
```

So, here, we just describe the whole problem, step by step:
- T is from a uniform distribution
- q is a deterministic function of T (and of r)
- k is a random variable whose distribution depends on q (and of N)
- k has been observed and is therefore 'clamped', or constrained, to be equal to this observed value k_obs.


Model variables
---------------------------------------------------
{:.subsubsection}


A very important point here: there are several types of variables in revbayes, which are defined using different assignment operators:

- constant variables, like $r$, or $N$ or $k_{obs}$.
- model variables: here, $T$, $q$ and $k$.
- MCMC variables: here, 'mymodel' (but the variables called 'moves', 'monitors' or 'mymcmc' below are also MCMC variables).

- constant variables are defined using the assignment operator '<-'
- the MCMC variables are defined using the assignment operator '='.
- as for model variables, they can be either stochastic ($T$ and $k$), or deterministic ($q$).
- stochastic model variables are defined using the stochastic operator '~'
- deterministic model variables are defined using the deterministic assignment operator ':='

- The model variables have the property that they remember how their value (for deterministic variables) or their probability distribution (for stochastic variables) depends on the other model variables. Thus, if other model variables change their value, then a model variable will 'know' it and will know how to recompute its value or its probability accordingly. This is very useful, because the whole idea of the MCMC is to 'move' the random variables (i.e. try small changes in their value) and accept or reject these moves depending on how this changes the global posterior distribution. Thus, each time we move a particular variable, all other variables will automatically know how they should update their value or probability, and the MCMC will be much more easily implemented.

- the value of a deterministic variable depends on the variables that are directly specified in its definition: here, $q$ depends on $T$ and $r$. Whenever the value of $T$ changes during the MCMC, this will change the value of $q$

- the probability of a stochastic variable depends on the variables that are directly specified in its definition: here, $k$ depends on $q$, and thus, if the value of $q$ changes during the MCMC (which will happen each time the value of $T$ will change), this will change the probability of observing $k = k_{obs}$.


The dependencies between the model variables can be visualized in terms of a graphical model. This point will be explained on the board.


Constructing the Monte Carlo inference
---------------------------------------------------
{:.subsubsection}


Consider the last line in the code show above:

```
mymodel = model(T)
```
With this instruction, the entire model will be captured by the RevBayes by searching for all stochastic and deterministic model variables that are interconnected, directly or indirectly, to T. Here, $q$ depends on $T$, and $k$ has a distribution that depends on $q$, so these three model variables, $T$, $q$, and $k$, will be captured in 'mymodel'.

Once the model is defined and conditioned (as done above), we should still specify how we want to 'move' the free random variables of the model. Here, we have clamped $k$, so it cannot move. As for $q$, it is not a random variable. On the other hand, $T$ is a random variable, and it is not clamped: it is in fact the variable that we want to estimate. So, this variable should be allowed to 'move' during the MCMC.
Here, we used a 'sliding' move:

```
moves[1] = mvSlide(T)
```
Note that 'moves' is a vector. Here, it has only one entry, but this is because our model has only one free random variable. In more complex models, we will have several types of moves, for different variables of the model (the tree, the rate of evolution, the nucleotide frequencies, the dN/dS, etc).

Next we want to monitor the MCMC. To do this, we make a vector of monitors (or monitoring streams):

```
monitors[1] = mnModel(filename = "apes_mcmc.log", printgen=10, separator = " ")
monitors[2] = mnScreen(printgen = 10, T)
```
The first monitor will write the sequence of values of T in a file called apes_mcmc.log. The second monitor will output a summary to the screen.

Finally, we construct an MCMC object, which gathers all of the modules: the model, the monitors, and the MCMC moves:

```
mymcmc = mcmc(mymodel, monitors, moves)
```

Now, we can run the model for 30 000 generations:
```
mymcmc.run(generations = 30000)
```
After the MCMC has run, a file names apes_mcmc.log has been created. This file contains all of the 3000 values of T visited during the MCMC (why only 3000?).

- Write this program, run it
- Draw the histogram of the values of $T$ visited during the MCMC
- Compare this histogram with the histogram obtained above with rejection sampling.


Monte Carlo inference using built-in phylogenetic tools
---------------------------------------------------
{:.subsubsection}


The model considered above is highly specialized and works only for two species.
We can re-implement the model considered above using the built-in functions and distributions proposed by RevBayes. This will give an introduction to the tools that we will then use in more complicated situations.

The idea is to proceed as follows:
- load the sequence alignment
- create model variables representing the tree and the age of the ancestor
- create model variables representing substitution process and its parameters
- create the substitution process and condition it on the sequence data
- implement MCMC moves
- run the MCMC

We now proceed step by step.
First, we load the data, from the nexus file called HC.nex, which contains only the sequences for Humans and Chimpanzees:
```
data <- readDiscreteCharacterData("HC.nex")
taxa <- data.taxa()
print(taxa)
```
Then, we specify the substitution rate (mutation rate per Myr), as above. Here, it is a constant:
```
# mutation rate per generation
mu <- 2e-8
# generation time
gentime <- 30
# mutation rate per Myr
r <- mu / gentime * 1e6
```
As above, we put a uniform prior on the age of the most recent common ancestor, between 0 and 20 Myr. And then, we create the tree. This is a time tree, which has only two branches, so we do not have to care about the tree topology.
```
T ~ dnUniform(0, 20)
tau ~ dnUniformTimeTree(T, taxa)
```
We then create the rate matrix for the Jukes-Cantor model
```
# JC substitution process
Q <- fnJC(4)
```
And then we invoke a substitution process along the tree tau, according to the Jukes Cantor substitution process Q, with an absolute substitution rate r, producing a DNA sequence alignment. Finally, we say that this process has led to the observed data; thus, we condition on the sequence alignment:
```
seq ~ dnPhyloCTMC(tree=tau, Q=Q, branchRates=r, type="DNA" )
seq.clamp(data)
```
As above, we pull out the model starting from one of the model variables (here, the age T):
```
mymodel = model(T)
```
And then we propose a scaling move on T, a series of monitors, and create the MCMC:
```
moves[1] = mvScale(T)

monitors[1] = mnModel(filename="jc.log", printgen=10, separator = TAB)
monitors[2] = mnScreen(printgen=100, T)

mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=30000)
```
- Write this program, run it
- Draw the histogram of the values of $T$ visited during the MCMC
- Compare this histogram with the histogram obtained above with rejection sampling.

