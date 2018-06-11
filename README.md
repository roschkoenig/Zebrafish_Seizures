# Analysing light-sheet imaging data of acute epileptic seizures in zebrafish 
_Code accompanying: Rosch et al (2018): Calcium-imaging and dynamic causal modelling reveal brain-wide changes in effective connectivity and synaptic dynamics during epileptic seizures. Accepted at PLoS Computational Biology_

This repository contains all code required to reproduce the dynamic causal modelling analysis of zebrafish light sheet imaging recordings during acutely induced epileptic seizures. Data together with the code are available through an OSF (open science framework) repository online. 

The code runs on [Matlab](https://uk.mathworks.com/products/matlab.html) (tested with 2016b), which unfortunately is not free, but often available through institutional subscriptions. The code also requires the following freely available software packages to run:

- [Statistical Parametric Mapping 12 (SPM12)](http://www.fil.ion.ucl.ac.uk/spm/) - this academic software implements the fundamental functions used for dynamic causal modelling (DCM) analysis
- [Color Brewer](https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab) - a great little package that produces perceptually balanced colour maps that are used for most of the plots shown here. A version of this toolbox is included in the repository. 

## Custom Routines included in this repository

The analysis is done by way of a number of custom routines, that in conjunction can be used to reproduce the findings from the paper. Below is a summary of what each of the routines does. Broadly, the code is trying to accomplish the following 5 objectives

1. Identify induced network-wide changes in neuronal activity at 'sensor space' - i.e. directly from the meausured signals
2. Use simulations to test whether DCM can resurrect neuronal parameters from calcium imaging signals
3. Use DCM and Bayesian model comparison to identify a parsimonious baseline network architecture
4. Use a hierarchical DCM to identify slow changes in neuronal parameters induced by the induced seizures 
5. Simulate the effects of the identified parameter changes to identify their effects on the network

### Visualise sensor space changes of neuronal dynamics using a sliding window
``zf_seizureexplore``

### Simulate dynamic causal modelling on the calcium dynamics 
``zf_calciumsim``

### Set up and invert baseline DCM 
``zf_dcm``

### Use Bayesian model reduction to make inference on model architecture at baseline
``zf_bmr``

### Set up (sliding window) files for DCM analysis 
``zf_slide`` and ``zf_slide_for_cluster`` 


### Use parametric empirical Bayes to make inference across time windows
``zf_peb``

### 

## Edits of corse SPM function called by the routines above
- ``zf_spm_dcm_csd_data``
- ``zf_spm_dcm_csd``
- ``zf_spm_dcm_fit``
- ``zf_spm_fs_csd``
- ``zf_spm_rand_power_law``
