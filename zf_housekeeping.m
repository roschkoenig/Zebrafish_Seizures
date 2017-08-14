function D = zf_housekeeping(sub)

% Housekeeping
%--------------------------------------------------------------------------
fs          = filesep;
if strcmp(computer, 'PCWIN64'), Fbase = 'C:\Users\rrosch\Dropbox\Research\Friston Lab\1608 Zebrafish';
else Fbase = '/Users/roschkoenig/Dropbox/Research/Friston Lab/1608 Zebrafish'; end

Fscripts    = [Fbase fs 'Scripts'];
Fanalysis   = [Fbase fs 'Matlab Files'];
Forig       = [Fbase fs 'Data'];
subs        = {'S1', 'S2', 'S3'};
lbl         = { 'RTect'; 'LTect'; 'RCrbl'; 'LCrbl'; ...
                'RRHbr'; 'LRHbr'; 'RCHbr'; 'LCHbr'; ... 
                'RRSpC'; 'LRSpC'};
            
addpath(genpath(Fscripts));

load([Forig fs 'single_plane_ROI_MEAN_TRACES']);
Z = ROI_MEAN_TRACES.data;
l = length(Z);

Fs          = 20;
fstps       = 200;
frq_ax      = linspace(1, Fs/2, fstps / 2);
win         = 60*Fs; 
stp         = 10*Fs;
windows     = 1:stp:l-win;


% Pack up for exporting
%==========================================================================
D.Fbase     = Fbase;
D.Fanalysis = Fanalysis;
D.Forig     = Forig;
D.Fscripts  = Fscripts;
D.lbl       = lbl;
D.subs      = subs;

D.Fs        = Fs;
D.frq_ax    = frq_ax;
D.stp       = stp;
D.win       = win;
D.endtime   = length(Z)/Fs/60;