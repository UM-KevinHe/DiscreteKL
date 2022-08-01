# DiscreteKL: Kullback-Leibler-Based Discrete Failure Time Models

This package is for the paper "Kullback-Leibler-Based Discrete Failure Time Models for Integration of
Published Prediction Models with New Time-To-Event Dataset"

## Purpose

Prediction models built based on a single data source may suffer from rare event rates, small sample sizes, and low signal-to-noise ratios. Incorporating published prediction models from large-scale studies is expected to improve the performance of prognosis prediction. Thus, we propose Kullback-Leibler-based (KL) discrete failure time models to integrate published prediction models (external models) with an individual-level time-to-event dataset (internal data) while accounting potential heterogeneity among different information sources. 

## Installation

```
#Install the package, need to install the devtools packages:
install.packages("devtools")
devtools::install_github("UM-KevinHe/DiscreteKL")

```

## Try it out with simulation examples

```
