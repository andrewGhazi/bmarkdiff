# bmarkdiff
Bayesian modelling of chipseq peaks

## Problem  

Running ChromHMM on tumor and normal samples, we get tens of thousands of genomic regions with identifiable marks. However, the number of samples in the cases and controls is much much smaller than the number of regions. So if we consider each region independently, there is not enough statistical power to identify regions that are differentially marked. 

## Motivation  

We noticed that there's a global correlation of counts in the two groups. We want this to serve as a background empirical prior. This provides estimate shrinkage, reducing the rate of false positives.

## Model

![Kruschke diagram](www/bmarkdiff_diagram.png)
