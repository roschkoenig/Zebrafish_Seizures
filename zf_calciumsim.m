%% Housekeeping
%==========================================================================
D           = zf_housekeeping(1);
fs          = filesep;

Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fanalysis   = D.Fanalysis;
Fdata       = [Fanalysis fs 'Data Files' fs 'Synth'];
Fdcm        = [Fanalysis fs 'DCM' fs 'Synth'];

spm('defaults', 'EEG');
addpath(Fscripts);
params = [];

if ~isempty(params)
    k = params.k;
    H = params.H;
    latency = params.latency;
else 
    k = 1/1000;   % Inverse time constant
    H = 1;        % Maximum height
    latency = 250;
end


%% Set up convolution kernel
%==========================================================================
clear cal ca

for ms = 1:10000
    ca(ms) = H * exp(-k*ms);
end

[val ind]   = max(ca);
uprise      = flip( (latency-1)^2 - [0:latency-1].^2 );
uprise      = uprise / max(uprise) * val;
cal = [uprise, ca];
ca  = cal(1:length(ca));

% Illustrate convolution kernel using spiking
% --------------------------------------------------------------------------
figure(1), clf

subplot(3,1,1), plot(ca);
rng('default')

Hz          = 1:50;
dt          = .001;
N           = 10000;

nsd         = [1./Hz]';
lfp         = real(zf_spm_rand_power_law(nsd,Hz,dt,N));

fc          = degtorad(360);
csrange     = linspace(0,fc,length(lfp));
mdt         = -cos(csrange) + 1;
intmt       = lfp .* mdt' + lfp;

subplot(3,1,2), plot(intmt);
con = conv(intmt, ca);
subplot(3,1,3), plot(con(1:N));


%% Model specification
%==========================================================================
range   = linspace(-1, 1, 6);
par     = 'H(1)';

sP = range;
rng('default')

% number of regions in simulated seizure activity and model specification
%--------------------------------------------------------------------------
Nc  = 1;                                   % number of channels
Ns  = 1;                                   % number of sources
options.spatial  = 'LFP';
options.model    = 'LFP';
options.analysis = 'CSD';
M.dipfit.model   = options.model;
M.dipfit.type    = options.spatial;
M.dipfit.Nc      = Nc;
M.dipfit.Ns      = Ns;
M.Hz             = 0.5:0.5:60;

% get associated priors
%--------------------------------------------------------------------------
A      = {0 0 0};
B      = {};
C      = 1;
pE     = spm_dcm_neural_priors(A,B,C,options.model);
pE     = spm_L_priors(M.dipfit,pE);
pE     = spm_ssr_priors(pE);
[x,f]  = spm_dcm_x_neural(pE,options.model);


% suppress channel noise (assuming many trials would be averaged)
%--------------------------------------------------------------------------
pE.a   = [ 0; 0];                  % log amplitude and f^(-a) exponent
pE.b   = [-2; 0];                  % log amplitude and f^(-a) exponent
pE.c   = [-2; 0];                  % log amplitude and f^(-a) exponent

% number of hidden states and endogenous inputs
%--------------------------------------------------------------------------
nx     = length(spm_vec(x));
nu     = size(pE.C,2);

% create LFP model
%==========================================================================
M.f    = f;
M.g    = 'spm_gx_erp';
M.x    = x;
M.n    = nx;
M.pE   = pE;
M.m    = nu;
M.l    = Nc;
M.u    = sparse(Ns,1);

P      = pE;

U.dt  = 1/500;
N     = 10/U.dt;
M.p   = 8;
M.dt  = U.dt;


clear nsd noise timsr convr
for t = 1:length(sP)
    
% Generate time series from predicted spectral densities
%--------------------------------------------------------------------------
switch par 
    case 'T(1)', P.T(1) = sP(t); 
    case 'T(2)', P.T(2) = sP(t); 
    case 'H(1)', P.H(1) = sP(t); 
    case 'H(2)', P.H(2) = sP(t); 
    case 'H(3)', P.H(3) = sP(t); 
    case 'H(4)', P.H(4) = sP(t); 
  	case 'H(5)', P.H(5) = sP(t);
end

nsd         = [1./M.Hz]';
noise(t,:)  = real(zf_spm_rand_power_law(nsd,M.Hz,U.dt,N));
noise(t,:)  = noise(t,:) / max(noise(t,:)) / 100;

psd         = spm_csd_mtf(P,M,U); 
timsr(t,:)  = real(zf_spm_rand_power_law(psd{1},M.Hz,U.dt,N));
timsr(t,:)  = timsr(t,:)/ max(timsr(t,:)) + noise(t,:);

tonvr       = conv(timsr(t,:), ca);
convr(t,:)  = tonvr(1:length(timsr(t,:)));

mar    = spm_mar(timsr(t,:)', M.p);
mar    = spm_mar_spectra(mar,M.Hz,1/U.dt);
csd    = mar.P;
clear mar

mar = spm_mar(convr(t,:)', M.p);
mar = spm_mar_spectra(mar, M.Hz, 1/U.dt);
ssd = mar.P;

end

timsr = timsr/max(max(timsr));
convr = convr/max(max(convr));

figure(2); clf;
for i = 1:size(timsr,1)    
subplot(length(sP), 1, i), 
    plot(timsr(i,:), 'k'); hold on
    plot(convr(i,:), 'r'); 
    
subplot(length(sP), 1, 1), 
    title('Simulated neuronal time series with calcium convolution');
end
    
% Make MEEG object with the time series
%==========================================================================
Fs = 1/U.dt;
ns = size(timsr,2);
nt = size(timsr,1);

ftdata.fsample  = Fs;
ftdata.label    = {'SynC'};
timaxis         = linspace(0, ns*U.dt, ns);

for s = 1:nt
    ftdata.trial{s} = timsr(s,:);
    ftdata.time{s}  = timaxis; 
    conds{s}        = ['LFP_' num2str(s)];
    
    ftdata.trial{nt + s}    = convr(s,:);
    ftdata.time{nt + s}     = timaxis;
    conds{nt+s}             = ['CAI_' num2str(s)];
end

cfg.resamplefs = 20;
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
cfg.feedback   = 'textbar';

ftdata = ft_resampledata(cfg, ftdata);

D = spm_eeg_ft2spm(ftdata, [Fdata fs 'syn.mat']);
D = type(D, 'single'); 
for c = 1:length(conds)
    D = conditions(D, c, conds{c});
end

S       = [];
S.task  = 'defaulteegsens';
S.D    	= D;
D       = spm_eeg_prep(S);
save(D);

Fs20Time = ftdata.time{1}*1000;

% Set up DCM structure and invert 
%==========================================================================
DCM = [];                           % create DCM struct 
DCM.options.analysis = 'CSD';       % cross-spectral density 
DCM.options.model    = 'LFP';      	% structure cannonical microcircuit (for now)
DCM.options.spatial  = 'LFP';        % virtual electrode input    
DCM.options.Tdcm     = [Fs20Time(1) Fs20Time(end)];   % 1-30k ms 
DCM.options.Fdcm     = [0.5 10];
DCM.options.D        = 1;         	% frequency bin, 1 = no downsampling

DCM.options.Nmodes   = 8;          	% cosine reduction components used 
DCM.options.han      = 0;         	% no hanning 


DCM.Sname           = {'SynC'};
DCM.Hz              = DCM.options.Fdcm(1):0.5:DCM.options.Fdcm(2);
DCM.xY.Dfile        = [Fdata fs 'syn.mat'];

DCM.A{1}        = 0;
DCM.A{2}        = 0;
DCM.A{3}        = 1;
DCM.B           = {};
DCM.C           = sparse(length(DCM.A{1}),0); 

% Define spatial inversion model (simple LFP)
%--------------------------------------------------------------------------
DCM.M.dipfit.Nm    = DCM.options.Nmodes;
DCM.M.dipfit.model = DCM.options.model;
DCM.M.dipfit.type  = DCM.options.spatial;
DCM.M.dipfit.Nc    = 1;
DCM.M.dipfit.Ns    = 1;

% Get priors
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);
DCM.pE   = pE;
DCM.pC   = pC;

% Run CAI
%--------------------------------------------------------------------------
for c = nt+1:length(conds)
    DCM.name = [Fdcm fs 'SynC_' conds{c}];
    DCM.options.trials  = [c];
    SYN{c} = zf_spm_dcm_fit(DCM);
end

for i = 1:(length(SYN)/2)
    CAI{i} = SYN{i+(length(SYN)/2)}{1};
end

% PEB to identify group mean
%==========================================================================
clear M;

% Set second level model parameters
%--------------------------------------------------------------------------
M.hE  = 0;
M.hC  = 1/16;
M.bE  = spm_vec(CAI{2}.M.pE);
M.bC  = diag(spm_vec(CAI{2}.M.pC));

M.X(:,1) = ones(1,length(CAI));
M.Q      = 'all';

% Run PEB
%--------------------------------------------------------------------------
CEB      = spm_dcm_peb(CAI', M, 'all');


% Rerun first level with new group means
%==========================================================================
% Equip first level DCM with new group mean priors
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);

vE              = spm_vec(pE);
vE(CEB.Pind)    = vE(CEB.Pind) + CEB.Ep(:,1);
pE              = spm_unvec(vE, pE);
DCM.M.pE   = pE;
DCM.M.pC   = pC;

clear CEB M CAI

% Run DCMs
%--------------------------------------------------------------------------
for c = nt+1:length(conds)
    DCM.name = [Fdcm fs 'SynC_' conds{c}];
    DCM.options.trials  = [c];
    SYN{c} = zf_spm_dcm_fit(DCM);
end

for i = 1:(length(SYN)/2)
    CAI{i} = SYN{i+(length(SYN)/2)}{1};
end


%% Run PEB for difference between conditions
%==========================================================================
clear M

%  Define model space in terms of parameters allowed to vary
%--------------------------------------------------------------------------
field = {'T(1)', 'T(2)', 'H(1)', 'H(2)', 'H(3)', 'H(4)', 'H(5)'};
M.X(:,1) = ones(1,length(CAI));
M.X(:,2) = [1:length(CAI)] - length(CAI)/2;

M.Q     = 'all';

for i = 1:length(field)
    [CEB{i} RCM{i}]  = spm_dcm_peb(CAI',M, field{i});
    CF(i,1) = CEB{i}.F;
end

[score winner] = max(CF);

BMC.CF  = CF;
BMC.CEB = CEB{winner};
BMC.CAI = CAI;

for r = 1:length(RCM{winner})
    Ep(r) = RCM{winner}{r}.Ep.H(1);
    Cps   = diag(RCM{winner}{r}.Cp);
    Cps   = spm_unvec(Cps, RCM{winner}{r}.Ep);
    Cp(r) = Cps.H(1);
end

figure
subplot(2,1,1), bar(BMC.CF - min(BMC.CF))
subplot(2,1,2), spm_plot_ci(Ep', Cp')


%% Plot convolution in frequency domain
%==========================================================================

lfp = timsr(1,:);

ftl = abs(fft(lfp)); fl = ceil(length(ftl)/2); 
ftl = ftl(1:fl);
ftl = smoothts(ftl, 'g', 50, 6);

fraxis = linspace(1/10, 250, fl);

dscal = cal(1:2:end);
tftc = abs(fft(dscal)); fl = ceil(length(tftc)/2);
tftc = tftc(1:fl);
ftc  = zeros(1,length(ftl));
ftc(1:length(tftc)) = tftc;
ftc     = ftc(1:length(ftl));

subplot(3,1,1), plot(fraxis, log(ftl));  xlim([0 20]);
subplot(3,1,2), plot(fraxis, log(ftc));  xlim([0 20]);
subplot(3,1,3), plot(fraxis, log(ftl .* ftc));    xlim([0 20]);

