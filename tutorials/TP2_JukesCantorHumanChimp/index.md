---
title: TP1 Implementing Bayesian inference with RevBayes
subtitle: Inferring the divergence time between Human and Chimpanzee
authors:  Nicolas Lartillot
level: 2
order: 0.2
index: true
prerequisites:
redirect: true
---

In this tutorial, we will first do simulations, and then conduct our first Bayesian analysis, using RevBayes. We have empirical data, in the form of a sequence alignment between two species (Human and Chimp). We can count the number of differences between the sequences of these two species. We also have some knowledge about the mutation rate in apes. Based on these data, we would like to infer the divergence time between Human and Chimp.

The data
===============
{:.section}


The file apes.nex contains an alignment with $N=878$ aligned nucleotide positions, between Human, Chimp and Gorilla. This alignment is a concatenation of 17 nuclear genes, for which only the 4-fold degenerate third codon positions have been kept. To visualize this alignment, you can open it with SeaView.
Out of these $N=878$ nucleotide positions, $k_{obs} =11$ are different between the human and the chimpanzee sequences.

We use only the 4-fold degenerate third codon positions because all mutations at these positions are synonymous, and therefore, we can assume that they are neutral (they are not under selection). As a result, the rate at which substitutions will accumulate over time is equal to the mutation rate.

Based on the whole-genome sequencing of families (trios), the mutation rate per generation in Human has been estimated at $u = 2.10^{-8}$. Similar experiments have been conducted in Chimpanzees, giving similar estimates. The generation time in Humans is around $\tau = 30$ years. It is a bit shorter in Chimpanzees, however, for simplicity, we will ignore this, and we will assume that the generation time has always been 30 years in great apes.

- Based on these estimates for $u$ and $\tau$, what is the mutation rate per million years? Let us call it $r$.

- Let us denote by $T$ the divergence time (measured in million years). Based on the observed number of differences between the two sequences (11 out of 878 positions), can you give a first rough estimate of this divergence time?


The Jukes-Cantor distance
===============
{:.section}

The rough estimate obtained in the last section does not take into account two complications:
- the problem of hidden substitutions: the fact that some positions may show the same nucleotide between the two species, and yet this is not the ancestral nucleotide (thus, there have been multiple substitutions along one or the two branches that have resulted in the same nucleotide being observed by chance). This problem is not so serious for closely related species such as Human and Chimp, but can be a more serious problem for more distantly related species.
- the fact that the total number of differences is small ($k_{obs} = 11$) and thus there is a large uncertainty about our estimate.

Here, we derive some mathematical identities that will give us a more accurate estimation procedure.

- First, given $r$ and $T$, we need to derive the probability that a given nucleotide position will be different between Human and Chimp.
Let us note this probability $q(r,T)$. It is a function of the rate $r$ and the divergence time $T$. Mathematically, it is given by:
$$$
q(r,T) = \frac{3/4} \, \left( 1 - e^{-2*4/3*r*T} \right)
$$$

- Second, given a sequence of length $N$, each of the $N$ positions has a probability of being different which is given by $q(r,T)$. As a result, the total number of differences for an alignment of $N=878$ positions, $k$, is random and has a binomial distribution:
$$
k \sim Binomial(N, q(r,T))
$$

We will now use these mathematical relations in two different ways. In a first step, we will simulate the evolution of genetic sequences, assuming a given value of the divergence time $T$ (in fact, we will only simulate the number of differences $k$, and not the entire DNA sequence). In a second step, we will rely on the observed number of differences $k_{obs} = 11$ to infer the divergence time $T$.

Simulating the evolutionary process
===============
{:.section}

In this section, we will use RevvBayes to write a little program that simulates the number of differences in the alignment, given the mutation parameters and the divergence time between the two species. This will give us a first practical introduction to RevBayes. We will do this directly in the command line environment offered by RevBayes. You can start revbayes by just typing:
```
rb
```
You should now be in the command line environment of RevBayes.

In this environment, we can first define the parameters of the model (using a syntax that looks very much like R):
```
# number of aligned nucleotide positions
N <- 878
# mutation rate per generation
u <- 2e-8
# generation time
tau <- 30
```
Based on these parameters, we can express the substitution rate per Myr:
```
# substitution rate per million years
r <- u / tau * 1e6
```

Then, for the purpose of simulation, we assume a certain value for the divergence time between human and chimp: let's say 6 Myr.
```
T <- 6
```
Now, we define a function which, given $r$ and $T$, returns the probability for a given position to be different between the two species (what we have noted $q(r,T)$ above):
```
q := Probability(( 1 - exp(-2*4/3*r*T)) * 3/4)
```
Here, there is an important thing to notice: we did not write 'q <- Probability...', but 'q := Probability'. What does that mean exactly? To understand this, you can make the following experiment: first, print the value of $q$:
```
print(q)
```
Now, change the value of $T$:
```
T <- 4
```
And finally, print again the value of $q$:
```
print(q)
```
Finally, change the value of $T$ back to 6, and print again $q$.
What do you see? Essentially, the variable $q$ 'keeps an eye' on the variable(s) on which it depends (here $T$). If those variables change their value, then $q$ also updates its value, based on the formula that was used for its definition.
Finally, we can proceed with the random variable $k$:
```
k ~ dnBinomial(q,N)
print(k)
```
Here, we use another symbol: '~', which is called the stochastic operator (as opposed to the deterministic assignment operators that we have seen above). This is just saying that $k$ is a random variable drawn from a binomial distribution of parameters $q$ and $N$.
Of note, if you try to print $k$ it several times (make the experiment), you will always get the same value (the one that was drawn when you first defined $k$). If you want to re-draw the variable, you have to do it explicitly:
```
k.redraw()
print(k)
```
In fact, you can use this to draw multiple times and see the distribution of values of $k$, using a for loop:
```
for (i in 1:10) { k.redraw(); print(k) }
```
Here again, you can try to change the value of $T$, say, set $T$ equal to 6, or to 20; and then redo 10 draws of $k$. What do you observe?


Bayesian inference and estimation of the divergence time
===============
{:.section}

(To avoid interference with the work done above, it can be useful to refresh the environment by closing the current RevBayes session, using the 'quit()' function, and then and open a new session.)

General overview:
To conduct Bayesian inference, we will use the same model as the one just explored, except that, now, we do not specify a value for $T$. Instead, we specify a prior distribution for $T$. This prior represents our state of knowledge about $T$, before we have seen the empirical data (the sequence alignment). In the present case, we will consider that we don't know anything a priori about this divergence time, except that it is between 0 and 20 Myr. Mathematically, we represent this by assuming a uniform distribution over $T$. In RevBayes, we can draw $T$ from this prior as follows:

```
T ~ dnUniform(0,20)
```
Thus, now, $T$ is a random variable (just like $k$). We just need to introduce this random $T$ in the model already considered above. So, we start by defining the fixed parameters (the constant) of the problem:

```
# number of aligned nucleotide positions
N <- 878

# mutation rate per generation
u <- 2e-8
# generation time
tau <- 30
# substitution rate per million years
r <- u / tau * 1e6
```
Next, we define our random $T$, from a uniform prior:
```
T ~ dnUniform(0,20)
```
Given $T$, we define, as before, the probability of a nucleotide difference $q$:
```
q := Probability(( 1 - exp(-2*4/3*r*T)) * 3/4)
```
and finally, the total number of differences $k$:
```
k ~ dnBinomial(N, q)
```


Now, we want to condition the model on the observed count (fix $k = 11$) and sample from the posterior distribution on $T$. To condition the model, we just  say that we have observed this random variable $k$, and it turns out that its value was 11 -- in other words, we 'clamp' it at this observed value:
```
k.clamp(11)
```

So, in summary, we can describe the whole problem, step by step:
- T is unknown and could be anything between 0 and 20 ($T$ is a random variable with a uniform distribution)
- q is a deterministic function of T (and of r)
- k is a random variable whose distribution depends on q (and of N)
- k has been observed and is therefore 'clamped', or constrained, to be equal to this observed value k_obs.
- what is the posterior distribution on T?


Constructing the Monte Carlo inference
---------------------------------------------------
{:.subsubsection}


Our model is now complete, and we capture it as a whole:
```
mymodel = model(T)
```
With this instruction, the entire model will be captured by RevBayes by searching for all stochastic and deterministic model variables that are interconnected, directly or indirectly, to T. Here, $q$ depends on $T$, and $k$ has a distribution that depends on $q$, so these three model variables, $T$, $q$, and $k$, will be captured in 'mymodel'.

Next, we need to specify how we want to 'move' the free random variables of the model. Here, we have clamped $k$, so it cannot move. As for $q$, it is not a random variable. On the other hand, $T$ is a random variable, and it is not clamped: it is in fact the variable that we want to estimate. So, this variable should be allowed to 'move' during the MCMC.
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
- Using Tracer, draw the histogram of the values of $T$ visited during the MCMC
- Compute the posterior mean, median and 95% credible interval for $T$

Model variables
---------------------------------------------------
{:.subsubsection}


A very important point here: there are several types of variables in revbayes, which are defined using different assignment operators:

- constant variables, like $r$, or $N$.
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

Of note, starting from now, it can be more convenient to write the series of commands, such as specified below, in a script, rather than in the command-line environment of RevBayes Then, you can run RevBayes directly on the script with the following command:
```
rb <scriptname>
```
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
