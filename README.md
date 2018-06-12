# Analysing light-sheet imaging data of acute epileptic seizures in zebrafish 
_Code accompanying: Rosch et al (2018): Calcium-imaging and dynamic causal modelling reveal brain-wide changes in effective connectivity and synaptic dynamics during epileptic seizures. Accepted at PLoS Computational Biology_

This repository contains all code required to reproduce the dynamic causal modelling analysis of zebrafish light sheet imaging recordings during acutely induced epileptic seizures. Data together with the code are available through an OSF (open science framework) repository online. 

The code runs on [Matlab](https://uk.mathworks.com/products/matlab.html) (tested with 2016b), which unfortunately is not free, but often available through institutional subscriptions. The code also requires the following freely available software packages to run:

- [Statistical Parametric Mapping 12 (SPM12)](http://www.fil.ion.ucl.ac.uk/spm/) - this academic software implements the fundamental functions used for dynamic causal modelling (DCM) analysis
- [Color Brewer](https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab) - a great little package that produces perceptually balanced colour maps that are used for most of the plots shown here. A version of this toolbox is included in the repository. 

## Custom Routines included in this repository

The analysis is done by way of a number of custom routines, that in conjunction can be used to reproduce the findings from the paper. Below is a summary of what each of the routines does. Broadly, the code is trying to accomplish the following 5 objectives

1. Identify induced [network-wide changes](#visualise-sensor-space-changes-of-neuronal-dynamics-using-a-sliding-window) in neuronal activity at 'sensor space' - i.e. directly from the meausured signals
2. Use simulations to test whether DCM can resurrect neuronal parameters from calcium imaging signals
3. Use [DCM](#set-up-and-invert-baseline-dcm) and [Bayesian model comparison](#use-bayesian-model-reduction-to-make-inference-on-model-architecture-at-baseline) to identify a parsimonious baseline network architecture
4. Use a [time-windowed](#set-up-sliding-window-files-for-dcm-analysis), hierarchical DCM to identify slow changes in neuronal parameters induced by the induced seizures 
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

These are then used as the basis to estimate the underlying neuronal parameter changes. For each of the (calcium-)time series, we will invert a single-node DCM for cross-specrtal densities. Using a Parametric empirical Bayesian appraoch, we will then identify the single parameter diesturbance that best explains the differences between these individual time traces. This will yield both a free energy approximation of the model evidence for any individual parameter to explain the observed effect (top panel of the figure below); and posterior estimates of the parameter change that best explains the transitions between the different time series. In this case, the correct parameter was identified from the PEB analysis (i.e. the one that was manipulated to simulate the original LFP traces).

![Bayesian model comparison and posterior parameter estimates](http://gdurl.com/voEb)

### Set up and invert baseline DCM 
*This will take a while to run*
``zf_dcm``

In the next step, we are using DCM to infer the model architecture that best explains the cross spectral density summaries of baseline network activity. For the DCM we are loading the baseline data for each animal, and setting up the model architecture and inversion parameters. 

```matlab 
DCM.options.analysis = 'CSD';       % cross-spectral density
DCM.options.model    = 'LFP';      	% three-population model
DCM.options.spatial  = 'LFP';       % virtual electrode input
```
This section defines the basic features of the model inversion - here we are analysing ongoing neuronal oscillatory activity using cross-spectral densities (CSD). We are assuming a basic structure that is currently described as the 'LFP' model - a three population summary of microcircuitry, consisting of a single main (excitatory) output population, and one inhibitory and one excitatory interneuronal population respectively. This is all analysed, treating the signal as direct neuronal recordings (i.e. a local field potential - LFP - spatial model).

The model space that is of interest here revolves around the 'A' matrix - i.e. the synaptic coupling between different neuronal sources. This is specified in `zf_modelspace`. Briefly, this designs a 2 by 2 by 6 model space defined by the presence or absence of forward/backward connectivity along the hierarchy; homologous connections between brain hemispheres; and additional extensive connections to and from a specified 'hub' region. At this stage, the full model is selected and used for inversion. 

### Use Bayesian model reduction to make inference on model architecture at baseline
``zf_bmr``

Based on the full model inversion for the baseline model performed at the preceding step, this function now fills in the remainder of the parameter estimates and approximated model evidence across the model space, using Bayesian model reduction, or BMR. This routine will set up a new DCM structure for each of the investigated models, then load the already inverted full DCM from the step above. Using the free energy approximation for the model evidence, this will then allow Bayesian model comparison, which the code will present as family-wise comparison as shown below. 

![Bayesian model reduction and comparison output](http://gdurl.com/Qc4d)

### Set up (sliding window) files for DCM analysis
*This takes a very long time to run, and will need setting up within your local infrastructure*
``zf_slide`` and ``zf_slide_for_cluster`` 

Based on the winning model architecture selected from the baseline model inversion performed in the previous two steps, we now invert DCMs for each individual time step separately - the idea here is that over short time period (i.e. 60s), a network of neural mass models is fitted to an assumed stationary oscillatory signal. By doing this with overlappying windows separatey by only short time steps, we will be able to track the slow parameter changes that underlie the transition between different neuronal states - such as here, the baseline and the induced seizure state. 

The `zf_slide` version of the code will specify the model for each time window and invert them in a single loop (over many days). The `zf_slide` version will specify the DCM structure needed for inversion, but not actually invert the DCMs. This will need to be done on a computing cluster with matlab installed, where a custom inversion function (`zf_spm_dcm_fit`) should be called. 

### Use parametric empirical Bayes to make inference across time windows
``zf_peb``

As we are interested in the changes of DCM parameters over time (i.e. between-DCM parameter changes), we can use parametric empirical Bayes (PEB) to specifically estimate between DCM effects that correspond to particular trajectories. For PEB we will specify a model space at the second (between-DCM) level that describes possible trajectories, and in the inversion we will identify the most parsimonious combination of effects of this second-level general linear model on DCM parameters to explain the data. 

Specifically in this instance the code will first load the DCMs that were inverted at the previous step. It will then define the types of between-DCM effects that we are looking for in terms of a second level design matrix in the variable `X`. For the purposes of this analysis, individual fish are treated as repeat measurements. We can visualise the design matrix by typing in `imagesc(X)` and should get the following outpout, where each row is a DCM / time window, and each column is an experimental effect or regressor. 

Similarly to the Bayesian model reduction performed above to identify a simpler baseline architecture, we can then use a similar approach define a set of models to compare at this (second) level - i.e. which of the model parameters are affected by the changes observed during seizure activity. We can then performe Bayesian model comparison between these reduced models, and plot the second level parameter estimates. 

### Simulate neuronal responses across parameter space
`zf_simulate`

Once we have identified a simple representation of consistent parameter changes across fish and time windows, we want to explore what effect those specific parameters have on the neuronal output. Because the DCM / PEB approach yields fully generative models that will produce predictions of neuronal oscillatory activty (in the shape of cross-spectral densities), we can simulate the effects of specific parameter changes.

This is illustrated in the paper using a single brain region as an example: For this we first use a low-dimensional representation of the data (through a principal component analysis of the windowed parameter estimates) and project each time window onto this low-dimensional representation. Having mapped out the space of the maximum variance, we can then simulate output spectra at each point, plotting the landscape of low-frequency power, but also high frequency power components not orignally included in the data (as this is a simulation). 


## Edits of corse SPM function called by the routines above
- ``zf_spm_dcm_csd_data`` - adapted to allow for smaller frequency bins (standard limit: 1Hz)
- ``zf_spm_dcm_csd`` - adapted to call `zf_spm_dcm_csd` 
- ``zf_spm_dcm_fit`` - adapted to call `zf_spm_dcm_csd`
- ``zf_spm_fs_csd`` - adapted to use log-scale of cross spectral densities
- ``zf_spm_rand_power_law`` - use additional scheme to generate random power law dynamics
