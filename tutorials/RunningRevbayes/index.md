---
title: Installing and using RevBayes with conda
subtitle:
authors:  Nicolas Lartillot, Bastien Boussau
level: 2
order: 0.8
prerequisites:
index: true
title-old: RB_Install
redirect: false
---


{% section Installing Revbayes %}

```
conda create --name revbayesenv
conda activate revbayesenv
conda config --add channels bioconda
conda config --add channels conda-forge
conda install  revbayes
conda install revbayes
conda install tracer
conda install figtree
conda deactivate
```



{% section Installing Revbayes on Pedago-NGS %}

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash miniconda.sh
conda init

conda config --add channels bioconda
conda config --add channels conda-forge
conda install revbayes
conda install tracer
conda install figtree

```


{% section Using Revbayes %}

```
conda activate revbayesenv
rb
conda deactivate
```

{% section Suggested folder organization %}

We like to create 3 folders: one for scripts, one for the data, and one for the output from the analyses.

```
mkdir scripts
mkdir analyses
mkdir data

```

To work on a script, one can run:

```
cd scripts
nano myScript.Rev

```
