% This code performs a sliding window analysis on the seizure data in order
% to visualise power distribution dynamics across the recording duration
%--------------------------------------------------------------------------

clear all

sliding  = 1;
plotting = 1;
saveplot = 0;

extrctBL = 0;
make_spm = 0;

% Housekeeping
%==========================================================================
fs          = filesep;
D           = zf_housekeeping;
Fbase       = D.Fbase;
Forig       = D.Forig;
Fanalysis   = D.Fanalysis;
lbl         = D.lbl;
subs        = D.subs;
scount = 0;


for sub = subs
scount = scount + 1
    
Fsub        = [Fanalysis fs 'Data Files' fs sub{:}];


switch sub{:}
    case 'S1', load([Forig fs 'single_plane_ROI_MEAN_TRACES']);
    case 'S2', load([Forig fs 'single_plane_s2_ROI_MEAN_TRACES']);
    case 'S3', load([Forig fs 'single_plane_s3_ROI_MEAN_TRACES']);
end

Z = ROI_MEAN_TRACES.data;

clear fullc fullt
Fs      = 20;
win     = 60 * Fs;
stp     = 10 * Fs;
l       = length(Z);
fstps   = win;
i       = 0;

if sliding
for s = 1:stp:l-win
    i               = i+1
    w               = Z(s:s+win, :);
    tfullt          = fft(w, fstps);
    fullt(i,:,:)   = abs(tfullt(1:floor(end/2),:));
    ft(i, :)        = mean(squeeze(fullt(i,:,:)), 1);
    
    fullc       = corr(w);
    halfc       = tril(fullc, -1);
    v           = halfc(find(halfc));
    co(i, :)    = v;
end

ccor{scount} = corr(co');
cpow{scount} = corr(abs(ft)');
end

% Plot spectral changes before and during seizure
%==========================================================================
% if plotting, 
cs = flip(cbrewer('div', 'RdGy', 100));
colormap(cs);
    
frq_axis = linspace(1, Fs, fstps / 2);
set(gcf, 'Color', 'w');

subplot(3,3,1:6)
    imagesc(log(squeeze(fullt(:,:,1)))'); hold on
    title('Log frequency-power over time', 'fontsize', 12, 'fontweight', 'bold');
    set(gca, 'YDir', 'normal');
    seg = floor(size(fullt,1)/6);

cp = cbrewer('qual', 'Paired', 10);

subplot(3,3,7)
    for i = 1:10
        plot(frq_axis, log(squeeze(mean(fullt(50:150,:,i),1))), 'Color', cp(i,:)); hold on
    end
    title('Pre PTZ average spectral distribution', 'fontsize', 12, 'fontweight', 'bold');
    ylim([4 14]);
    ylabel('Power');
    xlabel('Frequency');

subplot(3,3,8)
    for i = 1:10
        plot(frq_axis, log(squeeze(mean(fullt(300:400,:,i),1))), 'Color', cp(i,:)); hold on
    end
    title('Early seizure spectra', 'fontsize', 12, 'fontweight', 'bold');
    ylim([4 14]);
    ylabel('Power');
    xlabel('Frequency');
    
subplot(3,3,9)
    for i = 1:10
        plot(frq_axis, log(squeeze(mean(fullt(500:600,:,i),1))), 'Color', cp(i,:)); hold on
    end
    title('Late seizure spectra', 'fontsize', 12, 'fontweight', 'bold');
    ylim([4 14]);
    ylabel('Power');
    xlabel('Frequency');
    legend(lbl);
    
end

%% Generate SPM Files
%==========================================================================
if make_spm
tim_ax  = linspace(0, ((l / Fs)-1)/60, l);
i       = 0;
ftdata  = [];

% Sliding window to separate out time windows
%--------------------------------------------------------------------------
for s = 1:stp:l-win
    i 	= i+1;
    dw  = Z(s:s+win, :)';
    tw  = tim_ax(1:win);
    conds{i} = num2str(tim_ax(s));
    
    ftdata.trial{i}     = dw;
    ftdata.time{i}      = tw;
end

lbl = {'RTect'; 'LTect'; 'RCrbl'; 'LCrbl'; 'RRHbr'; 'LRHbr'; 'RCHbr'; 'LCHbr'; 'RRSpC'; 'LRSpC'};
ftdata.label = lbl;
ftdata.label = ftdata.label(:);

winstr  = num2str(floor(win/Fs));
stpstr  = num2str(floor(stp/Fs));

D = spm_eeg_ft2spm(ftdata, [Fsub fs 'Z_' winstr 'by' stpstr]);
D = type(D, 'single'); 
for c = 1:length(conds)
    D = conditions(D, c, conds{c});
end;

S = [];
S.task = 'defaulteegsens';
S.D = D;
D = spm_eeg_prep(S);

save(D);
end


%% Generate SPM Files
%==========================================================================
if extrctBL
    
l       = 4800;
tim_ax  = linspace(0, ((l / Fs)-1)/60, l);
ftdata  = [];


% Extract time window
%--------------------------------------------------------------------------
dw = Z(1:l,:)';
tw = tim_ax;

conds{1} = 'Baseline';
    
ftdata.trial{1}     = dw;
ftdata.time{1}      = tw;

lbl = {'RTect'; 'LTect'; 'RCrbl'; 'LCrbl'; 'RRHbr'; 'LRHbr'; 'RCHbr'; 'LCHbr'; 'RRSpC'; 'LRSpC'};
ftdata.label = lbl;
ftdata.label = ftdata.label(:);

D = spm_eeg_ft2spm(ftdata, [Fsub fs 'Baseline_data']);
D = type(D, 'single'); 
for c = 1:length(conds)
    D = conditions(D, c, conds{c});
end;

S = [];
S.task = 'defaulteegsens';
S.D = D;
D = spm_eeg_prep(S);

save(D);
end


end

cols = flip(cbrewer('div', 'Spectral', 100));
colormap(cols)

subplot(3,4,1), imagesc(cpow{1}, [0 1]); axis square
subplot(3,4,5), imagesc(cpow{2}, [0 1]); axis square
subplot(3,4,9), imagesc(cpow{3}, [0 1]); axis square

mcpow = zeros(size(cpow{1},1), size(cpow{1},2));
for cc = 1:length(cpow), mcpow = mcpow + cpow{cc}; end
mcpow = mcpow ./ length(cpow);

subplot(3,4,[2:4, 6:8, 10:12]);
imagesc(mcpow)
axis square

%% Plot Power dynamics and do statistics
%--------------------------------------------------------------------------
r1 = 1:240;    r2 = 241:400;   r3 = 401:800;

d{1} = [cpow{1}(1,r1), cpow{2}(1,r1), cpow{3}(1,r1)];
d{2} = [cpow{1}(1,r2), cpow{2}(1,r2), cpow{3}(1,r2)];
d{3} = [cpow{1}(1,r3), cpow{2}(1,r3), cpow{3}(1,r3)];

for dd = 1:length(d)
    c{dd} = ones(1,length(d{dd}));
end

zf_dotplot(d, {'Baseline', 'Seizure', 'Prolonged Seizure'}, 1, c, 0.1);
set(gcf, 'Color', 'w');
title('Power distribution correlations', 'fontsize', 12, 'fontweight', 'bold');
ylabel('Correlation with first time window');
