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

This code performs a sliding window analysis of the regionally averaged light sheet microscopy data. The code will look for the original data files (stored in 'Data' in the repository, and called something like this `single_plane_ROI_MEAN_TRACES.mat`). The code will then take 60s time windows in 10s steps to estimate time changing fourier spectra of the calcium signal and plot these (if specified). This function is also used to generate the SPM files (called MEEG objects) that will be required for the later analysis and are stored in the 'Matlab Files' folder. 

![Windowed spectral estimates](http://gdurl.com/2Dnf)

The windowed spectral estimates are also used to estimate a power dynamics correlation matrix (see discussion of the methods in [Rosch et al. 2018 Network Neuroscience, 2(1)](https://doi.org/10.1162/NETN_a_00026)) - for which the average across all fish is plotted as output, also seen in Figure 2 in the [bioRxiv preprint](https://doi.org/10.1101/160259). This matrix shows the correlation of each time window's channel-resolved frequency-power spectra with each other time window. The leading diagonal of the matrix is therefore always = 1 (i.e. fully correlated), and time periods where the overall output remains relatively static appear as blocks on the dynamics matrix. The function's output is shown below. 

![Power dynamics time by time correlation matrix](http://gdurl.com/J_Vn)

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
