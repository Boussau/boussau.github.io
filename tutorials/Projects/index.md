---
title: Projects using RevBayes
subtitle: 
authors:  Nicolas Lartillot, Bastien Boussau
level: 3
order: 1
prerequisites:
- All
index: true
title-old: RB_CTMC_Tutorial
redirect: false
---


{% section Coevolution between substitution rates and life-history traits %}

In mammals, there is substantial variation in substitution rate between lineages. This variation may be correlated with life-history traits: small mammals, with short generation times, tend to evolve faster. Here, the aim is measure the correlation between substitution rate and body mass or other life-history traits in mammals, by building a model of the joint evolutionary process of the substitution rate and life-history traits across branches, using a variant of the Brownian model used for modelling the auto-correlated molecular clock (Lartillot and Delsuc, 2012).

{% section A relaxed molecular clock informed by life-history and germ-line development in primates %}

In primates, the mutation rate per generation is mostly determined by the number of replications in the germline. The developmental process of the germline is relatively well characterized, and a model of its modulation, as a function of life-history (age of puberty, generation time, etc) has been proposed (Amster and Sella, 2016). Here, the aim is to use this model (and the empirically measured values reported in the article of Amster and Sella) to calibrate the molecular clock (and its variation) in simian primates, using information about sexual maturity and generation time in extant species.

{% section Estimating the strength of GC-biased gene conversion %}

The substitution patterns at neutral nucleotide positions in mammalian genomes is not just determined by the mutational process. Another force, GC biased gene conversion (gBGC) is also impacting substitution rates, by increasing the probability of fixation of G and C, compared to A and T. Using population genetics theory, it is possible to derive the nucleotide substitution rate matrix, as a function of the mutation rates and the intensity of gBGC (Lartillot, 2013). The aim is to develop this model and apply it to primate nucleotide data, so as to estimate the intensity of gBGC (and possibly, its variation between branches) using a nucleotide sequence alignment.

{% section Nucleotide compositional variation across species and its impact on phylogenetic reconstruction %}

The models that we have considered in the tutorials are homogeneous across branches. As a result, they predict the same nucleotide composition in all species. In practice, this is not the case. In some extreme situations, not accounting for compositional variation can result in phylogenetic reconstruction artifacts (typically, unrelated species with similar compositional biases artifactually cluster together). A simple but particularly striking example is shown in Foster, 2004. In this article, a model with branch specific nucleotide composition is also introduced. Here, the aim is to design a similar version of this model with branch-specific GC content, and to see if this improves phylogenetic reconstruction.

{% section Multi-gene phylogenetic reconstruction of the phylogeny of placental mammals. %}

The idea is to design a model for doing multi-gene phylogenetic reconstruction, in which all genes share the same species phylogeny but differ in their rate of evolution and in their GC content. The model can be used to reconstruct the phylogeny with greater accuracy (compared to single-gene analyses) but also, to estimate the variance in substitution rates and in GC content (or in other aspects of the substitution process) across genes. Two applications can be considered: to reconstruct the phylogeny of mammals (as in Murphy et al, 2001), or the phylogeny of simian primates.

{% section Molecular dating of Holocephali using fossil calibrations %}


