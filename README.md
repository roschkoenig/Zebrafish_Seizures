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
This code illustrates how parameters changes that correspond to underlying synaptic changes can be retrieved from both neurophysiological recordings, and calcium imaging data. For this, we are first simulating some LFP data, and illustrating a convolution with a kernel that resembles the temporal dynamics of a GCaMP6f probe (also used in the experimental paradigm). The results are shown below and will be given as output from this code - they show the LFP, the calcium kernel, and calcium imaging time series both in the temporal domain (left) and in the frequency domain (right). 

Temporal domain   | Frequency domain  
---               | ---          
Top to bottom: Calcium kernel; example LFP; convolved calcium time series | Top to bottom: LFP spetral density; Calcium kernel spectral density; resultant calcium time series spectral density 
![Calcium kernel in the time domain](https://gdurl.com/PVMc) | ![Calcium kernel in the frequency domain](https://gdurl.com/snWb)

The code will then generate a number of time series with a known set of neuronal parameters, where one of the parameters is varied across the different repeated simulations. 

![Several time series with varying parameters](https://gdurl.com/vUlE)

These are then used as the basis to estimate the underlying neuronal parameter changes - thus for each of the (calcium-)time series, we will invert a DCM for cross-spectral densities modelling a single node 

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
