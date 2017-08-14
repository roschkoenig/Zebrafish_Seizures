for sub = {'S1', 'S2', 'S3'};
%% Housekeeping
%--------------------------------------------------------------------------
clear DCM
fs          = filesep;
Fbase     	= '/Users/roschkoenig/Dropbox/Research/Friston Lab/1608 Zebrafish'; 
Fscripts    = [Fbase fs 'Scripts'];
Fanalysis   = [Fbase fs 'Matlab Files'];
Fdata       = [Fanalysis fs 'Data Files'];
Forig       = [Fbase fs 'Data'];
Fdcm        = [Fanalysis fs 'DCM' fs sub{:}];
datafile    = [Fdata fs sub{:} fs 'Baseline_data'];
spm('defaults', 'EEG');
addpath(Fscripts);

load([Forig fs 'single_plane_ROI_MEAN_TRACES']);
Z = ROI_MEAN_TRACES.data;
l = length(Z);

Fs          = 20;
fstps       = 200;
frq_ax      = linspace(1, Fs, fstps / 2);
tim_ax      = linspace(0, ((l / Fs)-1)/60, l);
win         = 60*Fs; 
stp         = 10*Fs;
windows     = 1:stp:l-win;
lbl = {'RTect'; 'LTect'; 'RCrbl'; 'LCrbl'; 'RRHbr'; 'LRHbr'; 'RCHbr'; 'LCHbr'; 'RRSpC'; 'LRSpC'};

%% Set up DCM structure and invert baselne
%==========================================================================
Amod = zf_modelspace;

DCM = [];                           % create DCM struct 
DCM.options.analysis = 'CSD';       % cross-spectral density 
DCM.options.model    = 'LFP';      	% structure cannonical microcircuit (for now)
DCM.options.spatial  = 'LFP';        % virtual electrode input    
DCM.options.Tdcm     = [tim_ax(1) tim_ax(end)];   % 1-30k ms 

DCM.options.Fdcm    = [frq_ax(1) frq_ax(end)];     	% frequency range  
DCM.options.Rft     = 5;           	% delay 
DCM.options.onset   = 64;         	% time delays between sources   
DCM.options.dur     = 16;       	% time dur  
DCM.options.D       = 0.2;         	% frequncy bin, 1 = no downsampling

DCM.options.Nmodes  = 8;          	% cosine reduction components used 
DCM.options.han     = 0;         	% no hanning 

DCM.options.lock    = 0;           	% lock the trial-specific effects  
DCM.options.multiC  = 0;            % multichannel effects  
DCM.options.location = 0;           % optmise location 
DCM.options.symmetry = 0;           % symmeterical dipoles

DCM.Sname           = lbl;
DCM.M.Hz            = frq_ax;
DCM.xY.Dfile        = datafile;

for a = 1:3
    DCM.A{a}    = Amod{end}.A{a};
end
DCM.B           = {};
DCM.C           = sparse(length(DCM.A{1}),0); 

DCM.M.dipfit.Nm    = DCM.options.Nmodes;
DCM.M.dipfit.model = DCM.options.model;
DCM.M.dipfit.type  = DCM.options.spatial;
DCM.M.dipfit.Nc    = 10;
DCM.M.dipfit.Ns    = 10;

[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);
DCM.pE   = pE;
DCM.pC   = pC;

DCM.name = [Fdcm fs 'Full_BLN'];
DCM.options.trials  = [1];
BLN = spm_dcm_fit(DCM);
end
% DCM.options.trials  = [64];
% PTZ = spm_dcm_fit(DCM);