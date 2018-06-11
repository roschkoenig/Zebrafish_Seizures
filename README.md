# Analysing light-sheet imaging data of acute epileptic seizures in zebrafish 
_Code accompanying: Rosch et al (2018): Calcium-imaging and dynamic causal modelling reveal brain-wide changes in effective connectivity and synaptic dynamics during epileptic seizures. Accepted at PLoS Computational Biology_

This repository contains all code required to reproduce the dynamic causal modelling analysis of zebrafish light sheet imaging recordings during acutely induced epileptic seizures. Data together with the code are available through an OSF (open science framework) repository online. 

The code runs on [Matlab](https://uk.mathworks.com/products/matlab.html) (tested with 2016b), which unfortunately is not free, but often available through institutional subscriptions. The code also requires the following freely available software packages to run:

- [Statistical Parametric Mapping 12 (SPM12)](http://www.fil.ion.ucl.ac.uk/spm/) - this academic software implements the fundamental functions used for dynamic causal modelling (DCM) analysis
- [Color Brewer](https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab) - a great little package that produces perceptually balanced colour maps that are used for most of the plots shown here. A version of this toolbox is included in the repo
