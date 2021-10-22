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



{% section Revbayes on Pedago-NGS %}

To work on Pedago-NGS, you first need to install conda. Once this has been done, installation of RevBayes can proceed as indicated further below.

To install conda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash miniconda.sh
conda init

```


{% section Installing Revbayes %}
We suggest that you install RevBayes in an *environment*. Environments are a nice way to encapsulate programs and settings without changing the parameters of your entire system. The environment we will create is called `revbayesenv`, and we will install all the software we need within this environment.

```
conda create --name revbayesenv
conda activate revbayesenv
conda config --add channels bioconda
conda config --add channels conda-forge
conda install revbayes
conda install tracer
conda install figtree
conda deactivate
```


{% section Using Revbayes %}

Assuming we have installed our programs in the environment `revbayesenv`, to use them we just need to activate it.

```
conda activate revbayesenv
rb
conda deactivate
```

{% section Suggested folder organization %}

We like to create 3 folders: one for scripts (`scripts`), one for the data (`data`), and one for the output from the analyses (`analyses`).

```
mkdir scripts
mkdir data
mkdir analyses

```

To work on a script, one can use a text editor like nano:

```
cd scripts
nano myScript.Rev

```
