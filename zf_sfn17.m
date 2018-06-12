% This code generates animations to be included in the dynamic poster
% presented at SfN 2017 in Washington DC
%--------------------------------------------------------------------------
% Housekeeping
%==========================================================================
clear all
fs          = filesep;
D           = zf_housekeeping;
Fbase       = D.Fbase;
Forig       = D.Forig;
Fanalysis   = D.Fanalysis;
lbl         = D.lbl;
subs        = D.subs;
scount = 1;

sub = subs(scount);

Fsub        = [Fanalysis fs 'Data Files' fs sub{:}];

switch sub{:}
    case 'S1', load([Forig fs 'single_plane_ROI_MEAN_TRACES']);
    case 'S2', load([Forig fs 'single_plane_s2_ROI_MEAN_TRACES']);
    case 'S3', load([Forig fs 'single_plane_s3_ROI_MEAN_TRACES']);
end

Z = ROI_MEAN_TRACES.data;

% Sliding window visualisation
%--------------------------------------------------------------------------
clear fullc fullt
Fs      = 20;
win     = 60 * Fs;
stp     = 10 * Fs;
l       = length(Z);
fstps   = win;
i       = 0;

for s = 1:stp:l-win
    i               = i+1
    w               = Z(s:s+win, :);
    tfullt          = fft(w, fstps);
    fullt(i,:,:)    = abs(tfullt(1:floor(end/2),:));
    ft(i, :)        = mean(squeeze(fullt(i,:,:)), 1);
    
    fullc       = corr(w);
    halfc       = tril(fullc, -1);
    v           = halfc(find(halfc));
    co(i, :)    = v;
end

ccor{scount} = corr(co');
cpow{scount} = corr(abs(ft)');

%%

h = figure;

set(gcf, 'position', [300 400 1300 600]);
imagecols   = flip(cbrewer('div', 'Spectral', 100));
timecols    = flip(cbrewer('div', 'Spectral', size(fullt,1)));
rawcols     = flip(cbrewer('qual', 'Paired', 10));
fax         = linspace(1/size(fullt,2), 10, size(fullt,2));

% Dynamic Matrix - STATIC
%--------------------------------------------------------------------------
subplot(2,4,[1 2 5 6]);
    colormap(imagecols); hold on
    imagesc(cpow{1});
    set(gcf, 'color', 'w');
    
%     corrbar = colorbar;
%     ylabel(corrbar, 'Correlation of windowed power spectra', 'fontsize', 12, 'rotation', 270);
%     set(corrbar, 'YTick', [-0.8 0.8], 'YTickLabel', [-0.8 0.8]);
    axis square

    title('Dynamic correlation matrix', 'fontsize', 16);
    set(gca, 'YTick', [1 fix(length(cpow{1})/5) length(cpow{1})], 'YTickLabel', {'0', 'PTZ', '150'})
    ylabel('Time [min]', 'fontsize', 14);
    
    set(gca, 'XTick', [1 fix(length(cpow{1})/5) length(cpow{1})], 'XTickLabel', {'0', 'PTZ', '150'})
    xlabel('Time [min]', 'fontsize', 14);
    
% Time series - STATIC
%--------------------------------------------------------------------------    
subplot(2,4, [3 4]);
    % Plots
    for z = 1:size(Z,2)
        plot([Z(:,z) - mean(Z(:,z))] + 1000*z, 'color', rawcols(z,:)); hold on
    end
       
    % Labels
    xlabel('time [3 min]', 'fontsize', 12);
    legend({'RTect'; 'LTect'; 'RCrbl'; 'LCrbl'; 'RRHbr'; 'LRHbr'; 'RCHbr'; 'LCHbr'; 'RRSpC'; 'LRSpC'});
    title('Regional average fluorescence');
   
    
for i = 1:size(fullt,1)
    s = 1:stp:l-win;
    
% Dynamic Matrix
%--------------------------------------------------------------------------
subplot(2,4,[1 2 5 6]);
    hold on
    if i > 1, oldy = newy; oldx = newx; end
    
    % Plotting Routines
    newy = plot([1 size(fullt,1)], [i i], 'w', 'linewidth', 3);
    newx = plot([i i], [1 size(fullt,1)], 'w', 'linewidth', 3);
    
    % Settings
    set(gca, 'Ydir', 'reverse');
    try delete(oldy); delete(oldx); end
    ylim([1 size(fullt,1)]);
    xlim([1 size(fullt,1)]);
 
% Time Series
%--------------------------------------------------------------------------    
subplot(2,4,[3 4]);
   
    % Settings
    xlim([s(i)-win s(i)+2*win]);  
    ylim([-Inf Inf]);
    set(gca, 'ycolor', [1 1 1]);
    set(gca, 'ytick', []);
    box off
    set(gca, 'xtick', []);
    
subplot(2, 4, [7 8]);
    % Plots
    p = plot(fax, log(squeeze(mean(fullt(i,:,:), 3))), 'color', cols(i,:)); hold on
    
    % Setting
    p.Color(4) = 0.2;
    xlim([-Inf Inf])
    ylim([4 12]);
    set(gca, 'ytick', []);
    box off
    
    % Labels
    title('Mean power spectra');
    xlabel('frequency [Hz]', 'fontsize', 12);
    ylabel('log power [a.u.]', 'fontsize', 12);
    
% Capture image as GIF
drawnow
frame   = getframe(h);
im      = frame2im(frame);
[imind cm] = rgb2ind(im, 256);

% Write to GIF File
if i == 1
    imwrite(imind, cm, 'test.gif', 'gif', 'Loopcount', Inf);
else
    imwrite(imind, cm, 'test.gif', 'gif', 'WriteMode', 'append');
end
end
