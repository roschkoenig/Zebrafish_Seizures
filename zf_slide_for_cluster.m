% Housekeeping
%--------------------------------------------------------------------------
clear DCM
fs          = filesep;
if strcmp('PCWIN64', computer); 
    Fbase = 'C:\Users\rrosch\Dropbox\Research\Friston Lab\1608 Zebrafish';
else
    Fbase = '/Users/roschkoenig/Dropbox/Research/Friston Lab/1608 Zebrafish';
end


for sub = {'S1', 'S2', 'S3'}
    
Fscripts    = [Fbase fs 'Scripts'];
Fanalysis   = [Fbase fs 'Matlab Files'];
Fdata       = [Fanalysis fs 'Data Files' fs sub{:}];
Forig       = [Fbase fs 'Data'];
Fdcm        = [Fanalysis fs 'DCM' fs sub{:}];
Fclust      = [Fanalysis fs 'Cluster Files' fs sub{:}];

Cbase       = ['/data/rrosch/zf_slide/'];
dataname    = 'Z_60by10.mat';
datafile    = [Fdata fs dataname];

spm('defaults', 'EEG');
addpath(Fscripts);

load([Forig fs 'single_plane_ROI_MEAN_TRACES']);
Z = ROI_MEAN_TRACES.data;
l = length(Z);
clear Z;

load(datafile);
Ntrials     = length(D.trials);
D.path      = Cbase;
save([Fclust fs dataname], 'D');
clear D;

Fs          = 20;
fstps       = 200;
frq_ax      = linspace(1, Fs, fstps / 2);
tim_ax      = linspace(0, ((l / Fs)-1)/60, l);
win         = 180*Fs; 
stp         = 90*Fs;
windows     = 1:stp:l-win;
lbl = {'RTect'; 'LTect'; 'RCrbl'; 'LCrbl'; 'RRHbr'; 'LRHbr'; 'RCHbr'; 'LCHbr'; 'RRSpC'; 'LRSpC'};

%% Set up DCM structure
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
DCM.xY.Dfile        = [Cbase dataname];
DCM.xY.Hz           = frq_ax;

for a = 1:3
    DCM.A{a}    = Amod{14}.A{a};     % Tectum Hub Model
end
DCM.B           = {};
DCM.C           = sparse(length(DCM.A{1}),0); 

DCM.M.dipfit.Nm    = DCM.options.Nmodes;
DCM.M.dipfit.model = DCM.options.model;
DCM.M.dipfit.type  = DCM.options.spatial;
DCM.M.dipfit.Nc    = 10;
DCM.M.dipfit.Ns    = 10;
DCM.xY.Ic          = DCM.M.dipfit.Nc;

[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);
DCM.pE   = pE;
DCM.pC   = pC;

for s = 1:Ntrials
    if s < 10, ns = ['00' num2str(s)]; 
    elseif s < 100, ns = ['0' num2str(s)];
    else ns = num2str(s); 
    end
    DCM.options.trials  = [s];
    save([Fclust fs 'DCM_for_Inversion' fs 'DCM_' ns '.mat'], 'DCM');
end
end