---
title: Two projects using RevBayes
subtitle:
authors:  Nicolas Lartillot, Bastien Boussau
level: 2
order: 0.8
prerequisites:
- EvolMolTP1
- EvolMolTP2
- EvolMolTP3
index: true
title-old: RB_Projects
redirect: false
---

{% section A relaxed molecular clock informed by [life-history and germ-line development in primates](PNAS2016Amster.pdf) %}

In primates, the mutation rate per generation is mostly determined by the number of replications in the germline. The developmental process of the germline is relatively well characterized, and a model of its modulation, as a function of life-history (age of puberty, generation time, etc) has been proposed (Amster and Sella, 2016). Here, the aim is to use this model (and the empirically measured values reported in the article of Amster and Sella) to calibrate the molecular clock (and its variation) in simian primates, using information about sexual maturity and generation time in extant species.

As a starting point, a script is provided, which implements a strict clock analysis using a multiple sequence alignment (only the 3d coding positions) for simians. This dating analysis is calibrated based on the mutation rate estimated in Humans. In addition, the program also models the evolution of the generation time along the tree, using a Gaussian model. You can run, and read, this script to understand its structure. Then, the idea is to use information about generation times along the tree to modulate the molecular clock along the branches, according to a (simplified version of) the mechanistic model suggested by Amster and Sella.

{% section Analysis of gene alignments from Gamma Proteobacteria %}

Here the project is to analyze alignments from several genes from gamma-Proteobacteria, some of which are free-living, and others are endosymbionts. We want to test whether selection is less efficient in endosymbiotic species, owing to their smaller effective population sizes.

Three single gene alignments for 30 bacteria are given, together with a pre-estimated tree and a file that contains the information about the status (free-living or endosymbiotic) for each of the 30 species. Finally, a script is given, gam30omega.rev, which implements a simple codon model that assumes the same dN/dS over all branches of the tree. Starting from this script, the idea is to implement a model that would assume a different dN/dS for the two types of bacteria.

