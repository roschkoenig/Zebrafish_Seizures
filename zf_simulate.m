%% Housekeeping
%==========================================================================
% Generic Housekeeping
%--------------------------------------------------------------------------
clear all
fs          = filesep;
D           = zf_housekeeping;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Forig       = D.Forig;
Fanalysis   = D.Fanalysis;
Fpeb        = [Fanalysis fs 'PEB'];

sub         = D.subs;
Fs          = D.Fs;
win         = D.win;
stp         = D.stp;
lbl         = D.lbl;
frq_ax      = D.frq_ax;
endtime     = D.endtime;

% Load outputs of PEB analysis
%--------------------------------------------------------------------------
Rlist  = cellstr(spm_select('FPList', Fpeb, '^R'));
for r = 1:length(Rlist)
    load(Rlist{r});
    RCM{r}  = rcm;
    clear rcm
end

load([Fpeb fs 'PEB']);
load([Fpeb fs 'PMA']);

%% Principal component analysis
%--------------------------------------------------------------------------
Np  = length(PEB.Pnames);
Nw  = length(PEB.Snames);
Ep  = [PEB.Ep(:,1:9) * PEB.M.X(:,1:9)']';

T = [];         H = [];
for p = 1:Np
    if ~isempty(regexp(PEB.Pnames{p}, ['T(.*[2],[1 2])'])),         T(p) = 1;    end
    if ~isempty(regexp(PEB.Pnames{p}, ['H(.*[2],[1 2 3 4 5])'])), 	H(p) = 1;    end  
end
T = find(T);    H = find(H);
colblock = flip(cbrewer('div', 'Spectral', Nw/3));
cols     = [colblock; colblock; colblock];

[Tcf Tsc]     = pca(Ep(:,T), 'Algorithm', 'eig');
[Hcf Hsc]     = pca(Ep(:,H), 'Algorithm', 'eig');

clear CA

CA.Tcf      = Tcf;          CA.Hcf      = Hcf;
CA.Tsc      = Tsc;          CA.Hsc      = Hsc;

scatter(Tsc(:,1), Hsc(:,1), 80, cols, 'filled')
xlabel('T component');
ylabel('H component');

%%
% 
% clear E Ep Er region regions
% tect    = [2];
% rspc    = [10];
% regions = {tect, rspc};
% H       = [];
% T       = [];
% 
% for rg = 1:length(regions)
% for r = 1:length(RCM)
% 
%     Na      = length(RCM{1}.Sname);
%     Ep      = RCM{r}.Ep;
%     pC   	= RCM{r}.pC; 
%     region  = regions{rg};
% 
%     % Find regional time constants
%     %----------------------------------------------------------------------
%     i = 1;
%     T = [T; Ep.T(region,:)];
%     
%     for t = 1:size(Ep.T,2)
%         E{i}    = Ep.T(region,t); 
%         i       = i + 1;
%     end
%     % 
%     % Find regional intrinsic connections
%     %----------------------------------------------------------------------
%     H = [H; Ep.H(region,:)];
%     for h = 1:size(Ep.H,2)
%         E{i}    = Ep.H(region,h);
%         i       = i + 1;
%     end
% 
%     % Restructure into single vector
%     %----------------------------------------------------------------------
%     Tp  = [];
%     for e = 1:length(E)
%         Tp  = [Tp; E{e}];
%     end
% 
%     Er{rg}(r,:)     = Tp;
%     
% end
% end
% 
% % Identify principal components
% %--------------------------------------------------------------------------
% [coeff score latent tsquared explained] = pca([Er{1}; Er{2}], 'Algorithm', 'eig');
% 
% [Tcf Tsc]   = pca(T, 'Algorithm', 'eig');
% [Hcf Hsc]   = pca(H, 'Algorithm', 'eig');
% 
% CA.T        = T;            CA.H        = H;
% CA.Tcf      = Tcf;          CA.Hcf      = Hcf;
% CA.Tsc      = Tsc;          CA.Hsc      = Hsc;
% 
% % Plot first two components
% %--------------------------------------------------------------------------
% r1  = 1:length(Er{1});
% r2  = r1 + length(Er{1});
% 
% % scatter(score(r1,1), score(r1,2), 40, 'filled', 'r'); hold on
% % scatter(score(r2,1), score(r2,2), 40, 'filled', 'b'); 
% 
% scatter(Tsc(r1,1), Hsc(r1,1), 40, 'filled', 'r'); hold on
% scatter(Tsc(r2,1), Hsc(r2,1), 40, 'filled', 'b'); 
% 
% 
% Np = length(PEB.Pnames);
% for rg = 1:length(regions)  
%     
% region  = regions{rg};
% rx      = '[';
% for rr  = 1:length(region)
%     if region(rr) > 9, sr = '0'; 
%     else               sr = num2str(region(rr)); end
%     rx = [rx sr ' '];
% end
% rx      = [rx ']'];
% 
% clear Te Ti He Hi
% 
% for p = 1:Np
%     if ~isempty(regexp(PEB.Pnames{p}, ['T.*' rx ',1)'])),          Te(p) = 1;    end
%     if ~isempty(regexp(PEB.Pnames{p}, ['T.*' rx ',2)'])),          Ti(p) = 1;    end  
%     if ~isempty(regexp(PEB.Pnames{p}, ['H.*' rx ',[1 2 3])'])),    He(p) = 1;    end  
%     if ~isempty(regexp(PEB.Pnames{p}, ['H.*' rx ',[4 5])'])),      Hi(p) = 1;    end 
% end
% 
% I(rg).Te    = find(Te);
% I(rg).Ti    = find(Ti);
% I(rg).He    = find(He);
% I(rg).Hi    = find(Hi);
% end
% 
% % clear acute prold
% % for rg = 1:length(region)
% %     idx             = [I(rg).Te I(rg).Ti I(rg).He, I(rg).Hi];
% %     acute(rg,:)     = [PEB.Ep(idx,1) + PEB.Ep(idx,2)]';
% %     prold(rg,:)     = [PEB.Ep(idx,3)]';
% %     
% %     P(rg).ac    = acute(rg,:) / coeff';
% %     P(rg).pl    = prold(rg,:) / coeff'; 
% % 
% %     sc = 20;
% %     plot([0 P(rg).ac(1)]*sc, [0 P(rg).ac(2)]*sc, 'linewidth', 4);
% %     plot([P(rg).ac(1) P(rg).pl(1)]*sc, [P(rg).ac(2) P(rg).pl(2)]*sc, 'linewidth', 4);
% %     axis square
% % end

%% Model single node output
%==========================================================================
frq_ax      = D.frq_ax;
BCM.A{1}    = 0;
BCM.A{2}    = 0;
BCM.A{3}    = 1;
BCM.B       = {};
BCM.C       = sparse(length(BCM.A{1}),0);

BCM.options.model    = 'LFP';

BCM.M.dipfit.Nm    = 8;
BCM.M.dipfit.model = BCM.options.model;
BCM.M.dipfit.type  = 'LFP';
BCM.M.dipfit.Nc    = 1;
BCM.M.dipfit.Ns    = 1;

frqsim      = .5:.5:50;
BCM.M.Hz  	= frqsim;
BCM.M.f     = 'spm_fx_lfp';
BCM.M.g     = 'spm_gx_erp';
BCM.M.l  	= BCM.M.dipfit.Nm;
BCM.M.x     = zeros(1,13);
BCM.M.u     = 0;
BCM.xU.X    = zeros(1,0);

[pE,pC]  = spm_dcm_neural_priors(BCM.A,BCM.B,BCM.C,BCM.options.model);
[pE,pC]  = spm_L_priors(BCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);

% Setup baseline 
%---------------------------------------------------------------------------
Ep  = PEB.Ep(:,7);
Np  = length(PEB.Pnames);

% clear T H
% for p = 1:Np
%     if ~isempty(regexp(PEB.Pnames{p}, ['T.*[2],[1 2])'])),          T(p) = 1;    end
%     if ~isempty(regexp(PEB.Pnames{p}, ['H.*[2],[1 2 3 4 5])'])), 	H(p) = 1;    end  
% end
% T   = find(T); 
% H   = find(H);

% Sp  = pE;
Sp.T    = Ep(T)';
Sp.H    = Ep(H)';

% 
% mnSc  = -50;   
% mxSc  = 30;    

steps = 100;
grad1  = linspace(-.5,.5,steps);
grad2  = linspace(-1,1,steps);
dim1  = grad1 .* CA.Tcf(:,1);
dim2  = grad2 .* CA.Hcf(:,1);

clear CSD POW
for d1 = 1:size(dim1,2)-1
for d2 = 1:size(dim2,2)
    pS      = Sp;
    ti      = [1 2];
    hi      = [1 2 3 4 5];
    pS.T    = pS.T + dim1(ti,d1);
    pS.H    = pS.H + dim2(hi,d2);
    
    csd         = spm_csd_mtf(pS, BCM.M, BCM.xU);
    CSD(d1,d2,:)= abs(csd{1}(:,1,1)); 
    POW(d1,d2)  = mean(CSD(d1,d2,1:8));
end
end

%%
for d1 = 1:size(dim1,2)-1
for d2 = 1:size(dim2,2)
    POW(d1,d2)  = mean(CSD(d1,d2,find(frqsim >= 40)));    
end
end
subplot(2,1,1)
contour(grad1(1:end-1), grad2, log(POW)', [-3 -2 -1 0 1]);
title('Log gamma power in parameter space');
colormap gray, axis square
colorbar;
set(gca, 'YDir', 'normal');

for d1 = 1:size(dim1,2)-1
for d2 = 1:size(dim2,2)
    POW(d1,d2)  = mean(CSD(d1,d2,find(frqsim <= 4)));    
end
end

subplot(2,1,2)
imagesc(grad1(1:end-1), grad2, log(POW)'); hold on
title('Log delta power in parameter space')
colormap gray, axis square
colorbar
set(gca, 'YDir', 'normal');

scatter(Tsc(:,1), Hsc(:,1), 80, cols, 'filled')
xlabel('T component');
ylabel('H component');


