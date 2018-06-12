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

sub         = D.subs;
Fs          = D.Fs;
win         = D.win;
stp         = D.stp;
lbl         = D.lbl;
frq_ax      = D.frq_ax;
endtime     = D.endtime;

i   = 0;
clear DCM

for s = 1:length(sub)

% Subject specific housekeeping
%--------------------------------------------------------------------------
Fdata       = [Fanalysis fs sub{s} 'Data Files'];
Fdcm        = [Fanalysis fs 'DCM' fs sub{s}];
Finv        = [Fanalysis fs 'Cluster Files' fs sub{s} fs 'Inverted DCMs'];
datafile    = [Fdata fs 'Z_60by10.mat'];


% Load DCMs 
%==========================================================================
files = cellstr(spm_select('FPList', Finv, 'DCM_*'));
for f = 5:5:892
    i = i + 1;
    TCM = load(files{f});
    DCM{i} = TCM.DCM{1};
    clear TCM;
end
end

DCM = DCM';

%% Set up PEB model
%==========================================================================
% Make second level model space
%==========================================================================
ldcm    = length(DCM)/3;
tim_ax  = linspace(0, endtime, ldcm);
ti      = find(tim_ax > 30);

% PTZ time curves
%--------------------------------------------------------------------------
clear ptz

k = 1/30;           % Inverse time constant
H = 1 * 1/0.37;     % Maximum height
i = 0;
ptz_t   = [];

for w = (tim_ax(tim_ax > 30))-30
    i = i+1;
    ptz(ti(i)) = H*k*w * exp(-k*w);
end

% Tonic effect
%--------------------------------------------------------------------------
tnc = zeros(1,ldcm);
tnc(ti) = ones(1,length(ti));

% Prolonged seizure effect
%--------------------------------------------------------------------------
plg = zeros(1,ldcm);
plg(ti) = linspace(0,1,length(ti));

% Direct cosine transforms
%--------------------------------------------------------------------------
dct = spm_dctmtx(1,4,linspace(0, 5, ldcm));
dct = dct(:,2:end)' ./ max(max(dct));

% Static effects
%--------------------------------------------------------------------------
stc = ones(1,ldcm);

% Empty filler
%--------------------------------------------------------------------------
emp = zeros(1,ldcm);

% Run PEB
%==========================================================================
clear X Xnames M PEB RCM BMA
X       = [ tnc, tnc, tnc; ... 
            ptz, ptz, ptz; ...
            plg, plg, plg; ...
            dct, dct, dct;
            stc, emp, emp; ...
            emp, stc, emp; ...
            emp, emp, stc]';
Xnames  = {'Tonic', 'Monophasic', 'Prolonged', 'DCT1', 'DCT2', 'DCT2', 'S1', 'S2', 'S3'};
M.X         = X;
M.Xnames    = Xnames;

imagesc(X);
[PEB RCM] = spm_dcm_peb(DCM, M, {'A', 'H', 'T'});

%% Create reduced model space at second level
%==========================================================================
clear T H

% Intrinsic effects
%==========================================================================
for t = 1:7, T{t} = zeros(10,2); H{t} = zeros(10,5); end
for f = 1:5
    idx     = [1:2]+(f-1)*2;
    T{f+1}(idx,:)   = DCM{1}.pC.T(idx,:);
    H{f+1}(idx,:)   = DCM{1}.pC.H(idx,:);    
end
T{7}      = DCM{1}.pC.T;        H{7}    = DCM{1}.pC.H;

% Extrinsic effects
%==========================================================================
Amod = zf_modelspace;
clear Hom Nei Hub hub nei extr

% Homologue lateral connectivity
%--------------------------------------------------------------------------
Hom{1} = Amod{1}.A;                             % none

% Neighbouring forward/backward connectivity
%--------------------------------------------------------------------------
empt    = zeros(length(Amod{1}.A{1}));
neimat  = zeros(length(Amod{1}.A{1}));

for m = 3:length(neimat)
    neimat(m, m-2) = 1;
    neimat(m-2, m) = 1;
end

Nei{1} = Amod{1}.A;                             % none
Nei{2} = Amod{1}.A;                             % forward
    Nei{2}{1} = empt + triu(neimat);                 
Nei{3} = Amod{1}.A;                             % backward
    Nei{3}{2} = empt + tril(neimat);
Nei{4} = Amod{1}.A;                             % forward and backward
    Nei{4}{1} = Nei{2}{1};
    Nei{4}{2} = Nei{3}{2};

% Optic Tectum hub forward/backward connectivity
%--------------------------------------------------------------------------
Hub{1} = Amod{1}.A;                             % none
Hub{2} = Amod{1}.A;                             % forward
    Hub{2}{1} = Amod{2}.A{1};
Hub{3} = Amod{1}.A;                             % backward
    Hub{3}{2} = Amod{2}.A{2};
Hub{4} = Amod{1}.A;                             % forward and backward
    Hub{4}{1} = Hub{2}{1};
    Hub{4}{2} = Hub{3}{2};

% Assemble into 4 * 4 = 32 model space
%--------------------------------------------------------------------------
count = 0;
h     = 1;
for n = 1:length(Nei)
for b = 1:length(Hub)
    count = count+1;
    ext{count}.A{1} = Nei{n}{1} + Hub{b}{1};
    ext{count}.A{2} = Nei{n}{2} + Hub{b}{2};
  	ext{count}.A{3} = Hom{h}{3};
    
    hub(count) = b;
    nei(count) = n;
end
end

% PEB BMC for Extrinsic Connections
%==========================================================================
% Assemble model matrix
%--------------------------------------------------------------------------
clear MOD

for tw = 1:length(DCM) 
    count = 0; 

    for i = length(T)
    for e = 1:length(ext)
        count = count+1;
        MOD{tw,count} = DCM{tw};

        MOD{tw,count}.M.pC.A     = ext{e}.A;
        MOD{tw,count}.M.pC.T     = T{i};
        MOD{tw,count}.M.pC.H     = H{i};
    end
    end
end

% Run PEB for Extrinsic Connections
%--------------------------------------------------------------------------
PMA     = spm_dcm_peb_bmc(PEB, MOD(1,:));

% Family Wise Comparison - Extrinsic Connections
%==========================================================================
figure
F = sum(PMA.F,2);
for ni = 1:length(Nei)
    n(ni) = sum(F(find(nei == ni)));
end
sn  = sort(n);
dFn = sn(end) - sn(end-1);

for hi = 1:length(Hub)
    h(hi) = sum(F(find(hub == hi)));
end
sh  = sort(h);
dFh = sh(end) - sh(end-1);

% Plot family-wise comparison for free energies 
%--------------------------------------------------------------------------
set(gcf, 'Color', 'w');

subplot(1,2,1)
bar(n - min(n))
    
    % Labels and Fonts
    title('Extrinsic coupling of hierarchy', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Nil', 'Forward', 'Backward', 'Forward/Backward'});
    xlabel(['dF = ' num2str(dFn,2)], 'FontWeight', 'bold');
    
    % Plotting Parameters
    axis square
    colormap gray

subplot(1,2,2)
bar(h - min(h))
    
    % Labels and Fonts
    title('Extrinsic coupling of hub region', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Nil', 'Forward', 'Backward', 'Forward/Backward'});
	xlabel(['dF = ' num2str(dFh, 4)], 'FontWeight', 'bold');
    
    % Plotting Parameters
    axis square
    colormap gray
    
    
%% PEB BMC for Intrinsic connections
%==========================================================================
% Assemble model matrix
%--------------------------------------------------------------------------
clear MOD t F

for tw = 1:length(DCM) 
    count = 0;
    for i = 1:length(T)
    for e = length(ext)
        count = count+1;
        
        MOD{tw,count} = DCM{tw};
        
        MOD{tw,count}.M.pC.A     = ext{e}.A;
        MOD{tw,count}.M.pC.T     = T{i};
        MOD{tw,count}.M.pC.H     = H{i};
        
        int(count) = i;
    end
    end
end

% Run PEB for Intrinsic Connections
%--------------------------------------------------------------------------
[PMA]   = spm_dcm_peb_bmc(PEB, MOD(1,:));

% Family Wise Comparison - Intrinsic Connections
%==========================================================================
clear t ti 
F   = sum(PMA.F,2);
st  = sort(F);
dFn = st(end) - st(end-1);

% Plot family-wise comparison for free energies 
%--------------------------------------------------------------------------
figure
bar(F - min(F))
set(gcf, 'Color', 'w');

% Labels and Fonts
title('Intrinsic Parameters', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', [1 2 3 4 5 6 7], 'XTickLabel', {'Nil', 'Tect', 'Crbl', 'RHBr', 'CHBr', 'RSC', 'all'});
xlabel(['dF = ' num2str(dFn,4)], 'FontWeight', 'bold');

% Plotting Parameters
axis square
colormap gray

%% Generate winning model average within model families
%==========================================================================
clear M
Xnames      = {'Tonic', 'Monophasic', 'Prolonged', 'DCT1', 'DCT2', 'DCT2', 'S1', 'S2', 'S3'};
M.X         = X;
M.Xnames    = Xnames;

Fhub        = { 'A{1}(1,3)', 'A{1}(1,4)', 'A{1}(1,5)', 'A{1}(1,6)', 'A{1}(1,7)', 'A{1}(1,8)', 'A{1}(1,9)', 'A{1}(1,10)', ...
                'A{1}(2,4)', 'A{1}(2,5)', 'A{1}(2,6)', 'A{1}(2,7)', 'A{1}(2,8)', 'A{1}(2,9)', 'A{1}(2,10)'};

[PEB, RCM]  = spm_dcm_peb(DCM, M, {'H', 'T', Fhub{:}});
[PMA]       = spm_dcm_peb_bmc(PEB);

%% Extract and plot Parameters
%==========================================================================
clear PEB PMA RCM

% Load outputs of PEB analysis
%--------------------------------------------------------------------------
Fpeb        = [Fanalysis fs 'PEB'];
load([Fpeb fs 'PEB']);
load([Fpeb fs 'PMA']);

Rlist  = cellstr(spm_select('FPList', Fpeb, '^R'));
for r = 1:length(Rlist)
    load(Rlist{r});
    RCM{r}  = rcm;
    clear rcm
end

%% Define time window to plot
%--------------------------------------------------------------------------
clear xEp xCp Te Ti He Hi Aa Ep Cp c d 
X         = PEB.M.X;
[val loc] = max(diff(PEB.M.X(:,2)));      % Maximum acute effect

G         = PEB;
for i = 1:2
    Np      = length(G.Pnames);
    Ep(:,i) = G.Ep( [1:Np] + (i-1)*Np ) * X(loc,i);
    Cps     = diag(G.Cp);
    Cp(:,i) = Cps( [1:Np] + (i-1)*Np );
end

xCp     = sum(Cp,2);
xEp     = sum(Ep,2);  xEp = exp(xEp);  

Te      = zeros(1,Np);  Ti  = zeros(1,Np);
He      = zeros(1,Np);  Hi 	= zeros(1,Np);
Aa      = zeros(1,Np);

for p = 1:Np
    if ~isempty(regexp(PEB.Pnames{p}, 'T.*,1)')), Te(p) = 1;    end
    if ~isempty(regexp(PEB.Pnames{p}, 'T.*,2)')), Ti(p) = 1;    end  
    if ~isempty(regexp(PEB.Pnames{p}, 'H.*,[1 2 3])')), He(p) = 1;  end  
    if ~isempty(regexp(PEB.Pnames{p}, 'H.*,[4,5])')), Hi(p) = 1;    end 
    if ~isempty(regexp(PEB.Pnames{p}, 'A.*')), Aa(p) = 1;    end 
end

Te = find(Te);  Ti = find(Ti);
He = find(He);  Hi = find(Hi);
Aa = find(Aa);

clear cols colblock
colblock = cbrewer('qual', 'Paired', 10);

d{1} = xEp(Te)';     lab{1} = 'Te';  c{1} = xCp(Te)'; cols{1} = [colblock];
d{2} = xEp(Ti)';     lab{2} = 'Ti';  c{2} = xCp(Ti)'; cols{2} = [colblock];
d{3} = xEp(He)';     lab{3} = 'He';  c{3} = xCp(He)'; cols{3} = [colblock; colblock; colblock];
d{4} = xEp(Hi)';     lab{4} = 'Hi';  c{4} = xCp(Hi)'; cols{4} = [colblock; colblock];
d{5} = xEp(Aa)';     lab{5} = 'A';   c{5} = xCp(Aa)'; cols{5} = [];

for a = 1:length(Aa)
    dgt         = str2double(PEB.Pnames{Aa(a)}(end-1));
    if dgt == 0, dgt = 10; end
    cols{5}     = [cols{5}; colblock(dgt,:)];
end

figure
zf_dotplot(d, lab, 1, c, cols);

