%% Housekeeping
%==========================================================================
fs          = filesep;
D           = zf_housekeeping;

Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fanalysis   = D.Fanalysis;
Forig       = D.Forig;
sub         = D.subs;

% Set up Bayesian Model Reduction across models of interest
%==========================================================================
% Set up BMR structure
%--------------------------------------------------------------------------
clear P
Amod = zf_modelspace;

% Collate all models in single matrix of models
%--------------------------------------------------------------------------
for s = 1:length(sub)
clear DCM 
Fdcm    = [Fanalysis fs 'DCM' fs sub{s}]; 
load([Fdcm fs 'Full_BLN.mat']);
FCM     = DCM;
DCM     = rmfield(DCM, 'Ep');
DCM     = rmfield(DCM, 'Cp');

for m = 1:(length(Amod)-1)
    for a = 1:3
        DCM.pE.A{a} = Amod{m}.A{a};
        DCM.pC.A{a} = DCM.pE.A{a} .* FCM.pC.A{a};
        DCM.M.pE    = DCM.pE;
        DCM.M.pC    = DCM.pC;
    end
    DCM.name    = Amod{m}.name;
    P{s,m}      = DCM;
end

% Add full model to model matrix
%--------------------------------------------------------------------------
P{s,length(Amod)} = FCM;

end

% Run Bayesian model reduction
%==========================================================================
[BLCM BLMC BLMA] = spm_dcm_bmr(P);

% Extract free energies
%--------------------------------------------------------------------------
clear F Fs
for b = 1:length(BLMC)
    F(b,:) = BLMC(b).F;
end

% Extract model names
%--------------------------------------------------------------------------
for m = 1:length(Amod)
    mlab{m} = Amod{m}.name;
end

%% Family wise comparison - draw figures
%==========================================================================
% Define families
%--------------------------------------------------------------------------
Flr = F(1:6) + F(7:12) + F(13:18) + F(19:24);
Fsr = F(1:6:24) + F(2:6:24) + F(3:6:24) + F(4:6:24) + F(5:6:24) + F(6:6:24);

figure(2)
set(gcf, 'color', 'w');
colormap gray

% Plot Hub connection BMC
%--------------------------------------------------------------------------
subplot(2,1,2), 
    % Plot
    bar(Flr - min(Flr));
    
    % Titles and legends
    title('Evidence for hub connections', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', 1:6);
    set(gca, 'XTickLabel', {'None', 'Tec', 'Cbl', 'RHb', 'CHb', 'RSC'});
    
    % Settings
%     axis square
    box off
    
% Plot Hub connection BMC
%--------------------------------------------------------------------------    
subplot(2,1,1), 

    % Plot
    bar(Fsr - min(Fsr));
    
    % Titles and Legends
    title('Evidence for short range connections', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', 1:4);
    set(gca, 'XTickLabel', {'None', 'Homol.', 'Neigh.', 'Both'})
    
    % Settings
%     axis square
    box off
    
