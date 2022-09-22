---
title: Projects using RevBayes
subtitle:
authors:  Nicolas Lartillot, Bastien Boussau
level: 2
order: 0.8
prerequisites:
- TP2_JukesCantorHumanChimp
- TP3_JC_PhylogenyReconstruction
- TP4_SubstitutionModels
- TP6_CodonModels
- TP7_RelaxedClock
index: true
title-old: RB_Projects
redirect: false
---


{% section Coevolution between [substitution rates and life-history traits](Evolution2012Lartillot-2.pdf) %}

In mammals, there is substantial variation in substitution rate between lineages. This variation may be correlated with life-history traits: small mammals, with short generation times, tend to evolve faster. Here, the aim is measure the correlation between substitution rate and body mass or other life-history traits in mammals, by building a model of the joint evolutionary process of the substitution rate and life-history traits across branches, using a variant of the Brownian model used for modelling the auto-correlated molecular clock ([Lartillot and Delsuc, 2012](Evolution2012Lartillot-2.pdf)).

{% section A relaxed molecular clock informed by [life-history and germ-line development in primates](PNAS2016Amster.pdf) %}

In primates, the mutation rate per generation is mostly determined by the number of replications in the germline. The developmental process of the germline is relatively well characterized, and a model of its modulation, as a function of life-history (age of puberty, generation time, etc) has been proposed (Amster and Sella, 2016). Here, the aim is to use this model (and the empirically measured values reported in the article of Amster and Sella) to calibrate the molecular clock (and its variation) in simian primates, using information about sexual maturity and generation time in extant species.

{% section Estimating the strength of [GC-biased gene conversion](MBE2013Lartillot-1.pdf) %}

The substitution patterns at neutral nucleotide positions in mammalian genomes is not just determined by the mutational process. Another force, GC biased gene conversion (gBGC) is also impacting substitution rates, by increasing the probability of fixation of G and C, compared to A and T. Using population genetics theory, it is possible to derive the nucleotide substitution rate matrix, as a function of the mutation rates and the intensity of gBGC (Lartillot, 2013). The aim is to develop this model and apply it to primate nucleotide data, so as to estimate the intensity of gBGC (and possibly, its variation between branches) using a nucleotide sequence alignment.

{% section Nucleotide [compositional variation across species and its impact on phylogenetic reconstruction](SystBiol2004Foster-1.pdf) %}
The models that we have considered in the tutorials are homogeneous across branches. As a result, they predict the same nucleotide composition in all species. In practice, this is not the case. In some extreme situations, not accounting for compositional variation can result in phylogenetic reconstruction artifacts (typically, unrelated species with similar compositional biases artifactually cluster together). A simple but particularly striking example is shown in [Foster, 2004](SystBiol2004Foster-1.pdf). In this article, a model with branch specific nucleotide composition is  introduced to improve phylogenetic reconstruction in such cases. Since then, [Heaps et al.](Heaps_2014.pdf) have proposed improved models. Here, the aim is to design versions of these models with branch-specific equilibrium frequencies, and to see if they improve phylogenetic reconstruction.


{% section [Multi-gene phylogenetic reconstruction of the phylogeny of birds](Jarvis_2014.pdf) %}
Designing a model for doing multi-gene phylogenetic reconstruction. Genes may share the same species phylogeny but may differ in their rate of evolution and in their GC content. Can be used to reconstruct the phylogeny or to estimate the variance in substitution rates and in GC content across genes. The data is available from [the article](Jarvis_2014.pdf)

{% section [Correlation between GC composition of ribosomal RNA and growth temperature in Archaea](Groussin_2011.pdf). %}

The idea is to model the correlated evolution of rRNA GC content and growth temperature across Archaea, and use this model to estimate the correlation and to infer ancestral temperatures along the phylogeny.

{% section Are patterns of absence/presence of genes across genomes informative about the phylogeny? %}

The idea is to model the process of gain and loss of genes across a phylogeny, and to apply this model to data of absence/presence of genes across metazoans. See article of [Pisani et al., (2015)](Pisani_2015.pdf) and [Ryan et al (2013)](Ryan_2013.pdf) for two analyses giving different results on the same dataset.


{% section [Convergent evolution toward subterranean lifestyle in isopods, and its consequences on the rate of genome evolution](Saclier_2018.pdf). %}

Saclier et al have analysed the evolutionary patterns in a group of isopod species, in which there has been a large number of independent transitions from surface to underground lifestyle -- giving an opportunity for modeling and investigating the impact of these transitions on genomic sequence evolution.


{% section Exploring the feasibility of [importance sampling](https://en.wikipedia.org/wiki/Importance_sampling) for phylogenetic inference %}

Bayesian phylogenetic inference classically relies on MCMC, a computationally intensive algorithm. This is a problem as data sets have been increasing in size, resulting in increased computational and environmental footprints. In this project, you will investigate an alternative algorithm, [importance sampling](https://en.wikipedia.org/wiki/Importance_sampling), as follows. Firstly, you will perform MCMC inference on a single gene alignment, and on a big alignment. Secondly, you will subsample a limited number of samples from the posterior distribution obtained on the small alignment. Thirdly, you will evaluate the posterior probability of these samples according to the big alignment. Fourthly, you will reweigh these samples according to the ratio of their posterior probabilities, on the large data set *vs* the small one. As a result, on the big dataset, you will then have a posterior distribution obtained through importance sampling, and one obtained using MCMC. On the big dataset, does the faster importance sampling approach produce estimates that are similar to the MCMC approach?


{% section Analysis of the SEMG2 gene in Primates %}

Here the project is to analyze an alignment of the SEMG2 gene from primate species that differ in their mating systems.


{% section Analysis of mitochondrial protein evolution in Daphnia %}

Here the project is to analyze a concatenated alignment of 15 DNA sequences coding for proteins from 29 strains of *Daphnia pulex*, some of which reproduce sexually (named S1 to S14), and others, asexually (named A1 to A14). Sexual reproduction is assumed to be the ancestral condition.
