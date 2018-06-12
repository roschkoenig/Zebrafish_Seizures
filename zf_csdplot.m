% ZF CSD Plot
%==========================================================================
% This function loads inverted DCMs of a specified subject and plots
% example CSDs to illustrate the goodness of fit of the DCM inversion

% Manual Definitions
%--------------------------------------------------------------------------
s           = 1;

% Housekeeping
%--------------------------------------------------------------------------
clear DCM
fs          = filesep;
D           = zf_housekeeping;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Forig       = D.Forig;
Fanalysis   = D.Fanalysis;

sub         = D.subs;
Fs          = D.Fs;
win         = D.win;
stp         = D.stp;
lbl         = D.lbl;
frq_ax      = D.frq_ax;

Fdata       = [Fanalysis fs sub{s} 'Data Files'];
Fdcm        = [Fanalysis fs 'DCM' fs sub{s}];
Finv        = [Fanalysis fs 'Cluster Files' fs sub{s} fs 'Inverted DCMs'];
datafile    = [Fdata fs 'Z_60by10.mat'];

load([Forig fs 'single_plane_ROI_MEAN_TRACES']);
Z = ROI_MEAN_TRACES.data;
l = length(Z);

% Load DCMs and extract obs and pred CSDs
%==========================================================================
files = cellstr(spm_select('FPList', Finv, 'DCM_*'));
for f = 1:length(files)
    TCM = load(files{f});
    DCM(f) = TCM.DCM{1};
    clear TCM;
end

for f = 1:length(files)
    Hc(f,:,:,:)     = DCM(f).Hc{1};
    y(f,:,:,:)      = DCM(f).xY.y{1};
end

% Plot example modes of model prediction
%==========================================================================
clear aHc tHc mHc ay ty my

% Extract diagonal modes
%--------------------------------------------------------------------------
aHc = abs(Hc);
for i = 1:8, tHc(:,:,i) = squeeze(aHc(:,:,i,i)); end
mHc = squeeze(mean(tHc,3));
    
ay = abs(y);
for i = 1:8, ty(:,:,i) = squeeze(ay(:,:,i,i)); end
my = squeeze(mean(ty,3));

%% Set up plotting variables
%--------------------------------------------------------------------------
cls = flip(cbrewer('seq', 'YlGnBu', 100));
colormap(cls);
set(gcf, 'Color', 'w');
tim_ax = linspace(1, length(Z)/Fs/60, length(Z));
stp_ax = linspace(1, tim_ax(end), size(Hc,1));

md      = [1 5];
plrange = [-4.5 2.5];

% Plot Mode 5
%--------------------------------------------------------------------------
subplot(6,1,1), imagesc(stp_ax, frq_ax, log(ty(:,:,md(1)))', plrange)
    set(gca, 'xtick', []);
    set(gca, 'ydir', 'normal');
    box off
    set(gca, 'tickdir', 'out')
    title('Example Mode 1: Observed', 'fontsize', 12, 'fontweight', 'bold')

subplot(6,1,2), imagesc(stp_ax, frq_ax, log(tHc(:,:,md(1)))', plrange);
    set(gca, 'xtick', []);
    set(gca, 'ydir', 'normal');
    box off
    set(gca, 'tickdir', 'out')   
    title('Example Mode 1: Predicted', 'fontsize', 12, 'fontweight', 'bold')

% Plot Time Series
%--------------------------------------------------------------------------
subplot(6,1,3), 
    plot((Z(:,6)), 'Color', [0.3 0.3 0.3]);
    axis off;

% Plot Mode 1
%--------------------------------------------------------------------------
subplot(6,1,4), imagesc(stp_ax, frq_ax, log(ty(:,:,md(2)))', plrange)
    set(gca, 'ydir', 'normal');
    set(gca, 'xtick', []);
    box off
    set(gca, 'tickdir', 'out')
    title('Example Mode 2: Observed', 'fontsize', 12, 'fontweight', 'bold')
    
subplot(6,1,5), imagesc(stp_ax, frq_ax, log(tHc(:,:,md(2)))', plrange);
    set(gca, 'ydir', 'normal');
    box off
    set(gca, 'tickdir', 'out')
    xlabel('time [min]');
    title('Example Mode 2: Predicted', 'fontsize', 12, 'fontweight', 'bold')

prd = squeeze(mean(log(aHc(:,:,:,:))));   prd = prd(:);
obs = squeeze(mean(log(ay(:,:,:,:))));    obs = obs(:);

t = table(prd, obs, 'variablenames', {'prd', 'obs'});
lm = fitlm(t, 'linear')

%% Code Grave

% % Plot example CSDs pre and during seizure activity
% %--------------------------------------------------------------------------
% cols = parula(length(DCM));  %cbrewer('div', 'Spectral', length(DCM));
% pr   = 1:100;
% sz   = 351:450;
% for md = [1 5]
% figure(md)
% set(gcf, 'Color', 'w');
% 
% for i = pr
%     subplot(2,2,1)
%     plot(frq_ax, squeeze(log(abs(Hc(i,:,md,md))))', 'Color', [0.8 0.8 0.9]); hold on;
%     subplot(2,2,2)
%     plot(frq_ax, squeeze(log(abs(y(i,:,md,md))))', 'Color', [0.7 0.7 0.8]); hold on;    
% end
% subplot(2,2,1), 
%     plot(frq_ax, mean(squeeze(log(abs(Hc(pr,:,md,md)))),1), 'Color', [0.2 0.2 0.4], 'linewidth', 2)
%     ylim([-8 4]);
%     ylabel('Log Power', 'fontweight', 'bold');
%     xlabel('Frequency', 'fontweight', 'bold');
%     title('Baseline: Predicted', 'fontsize', 12, 'fontweight', 'bold');
% subplot(2,2,2), 
%     plot(frq_ax, mean(squeeze(log(abs(y(pr,:,md,md)))),1), 'Color', [0.1 0.1 0.3], 'linewidth', 2)
%     ylabel('Log Power', 'fontweight', 'bold');
%     xlabel('Frequency', 'fontweight', 'bold');
%     title('Baseline: Observed', 'fontsize', 12, 'fontweight', 'bold');
% 	ylim([-8 4]);
% 
% 
% for i = sz
%     subplot(2,2,3)
%     plot(frq_ax, squeeze(log(abs(Hc(i,:,md,md))))', 'Color', [0.8 0.9 0.8]); hold on;
%     subplot(2,2,4)
% 	plot(frq_ax, squeeze(log(abs(y(i,:,md,md))))', 'Color', [0.7 0.8 0.7]); hold on; 
% end
% subplot(2,2,3), 
%     plot(frq_ax, mean(squeeze(log(abs(Hc(sz,:,md,md)))),1), 'Color', [0.2 0.4 0.2], 'linewidth', 2)
%     title('Seizure: Predicted', 'fontsize', 12, 'fontweight', 'bold');
%     ylabel('Log Power', 'fontweight', 'bold');
%     xlabel('Frequency', 'fontweight', 'bold');
% 	ylim([-8 4]);
% subplot(2,2,4),    
%     plot(frq_ax, mean(squeeze(log(abs(y(sz,:,md,md)))),1), 'Color', [0.1 0.3 0.1], 'linewidth', 2)
%     title('Seizure: Observed', 'fontsize', 12, 'fontweight', 'bold');
% 	ylabel('Log Power', 'fontweight', 'bold');
%     xlabel('Frequency', 'fontweight', 'bold');
% 	ylim([-8 4]);  
% end
