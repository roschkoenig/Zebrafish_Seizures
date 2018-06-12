function zf_dotplot(d, lab, sig, cv, varargin)
%==========================================================================
% This function takes data in a cell array and plots each dataset within a
% cell as individual dotplots. It then performs a 2 sided ttest on all
% variations and plots means and standard deviations, as well as the
% significance of the mean differences

% Unpack function parameters
%--------------------------------------------------------------------------
if nargin > 4, 
    colm    = varargin{1}; 
    usecol  = 1;
else
    usecol  = 0;
end

scl = 1;

% Identify range of covariances
%--------------------------------------------------------------------------
for cc = 1:length(cv)
    mx(cc) = max(cv{cc});
    mn(cc) = min(cv{cc});
end
cvrange = [min(mn) max(mx)];

% Loop through parameter groups to be plotted
%==========================================================================
for dd = 1:length(d)
clear a cols colsort TCM 

% Define offset jitter
%--------------------------------------------------------------------------
o       = ones(1,length(d{dd})) * dd + randn(1,length(d{dd}))/20;

% Set colours
%--------------------------------------------------------------------------
if usecol,      cols    = colm{dd};
else
    gs  = .5;
    for di = 1:length(d{dd}), cols(di,:) = [gs gs gs]; end
end

% Define scaling factor to represent covariance of the estimates
%--------------------------------------------------------------------------
allcvs = [cv{:}];
for di = 1:length(d{dd})
    ci      = ceil((cv{dd}(di) - min(cv{dd}))/(range(cv{dd})) * 70)+1;
    a(di) 	= 101 - ci;
end

% Define bar length of individual confidence intervals
%--------------------------------------------------------------------------
% for di = 1:length(d{dd})
%     conf    = spm_invNcdf(1-.1); 
%     c(di)   = conf*sqrt(cv{dd}(di));
% end

% Plot
%==========================================================================
% for di = 1:length(d{dd})
%     up = exp(log(d{dd}(di)) + c(di));
%     lo = exp(log(d{dd}(di)) - c(di));
%     plot([o(di) o(di)], [lo up], 'color', cols(di,:)); hold on
% end
scatter(o, d{dd}, a, cols, 'filled'); hold on;
xlim([0, length(d)+1]);

% Estimate Bayesian parameter averages if indicated
%--------------------------------------------------------------------------
if sig == 2
for di = 1:length(d{dd})
    TCM{di}.Ep      = d{dd}(di);
    TCM{di}.Cp      = cv{dd}(di);
    TCM{di}.M.pE    = 0;
    if strcmp(lab{dd}, 'A'),    TCM{di}.M.pC    = .5;
    else                        TCM{di}.M.pC    = .065;     end
end
BPA = spm_dcm_bpa(TCM);
end

% Define confidence intervals or ranges to be plotted
%--------------------------------------------------------------------------
if sig < 2
    m   = median(d{dd});
    upp = quantile(d{dd},.75);
    low = quantile(d{dd},.25);

elseif sig == 2
    m       = BPA.Ep;
    civl    = spm_invNcdf(1 - 0.05);
    civl    = civl * sqrt(BPA.Cp);
    upp     = BPA.Ep + civl;
    low     = BPA.Ep - civl;
    Ees(dd)     = BPA.Ep;
    Ccs(dd)     = BPA.Cp;
end

% Plot mean estimates and ranges
%--------------------------------------------------------------------------
plot([dd dd], [low upp], 'Color', [0.7 0.3 0.3], 'Linewidth', 2);
plot([dd-0.25 dd+0.25], [m m], 'k', 'Linewidth', 2);

end

% Define plotting parameters
%--------------------------------------------------------------------------
set(gca, 'XTick', 1:length(d));
set(gca, 'XTickLabel', lab, 'FontWeight', 'bold');

l = combnk(1:length(d), 2);
i = 0;
t = max([d{:}]); st = 0.1 * t;

% Plot significance tests
%--------------------------------------------------------------------------
if sig == 1
for ll = 1:size(l,1)
    p  = ranksum(d{l(ll,1)}, d{l(ll,2)});
    if p < 0.05
       plot([l(ll,1) l(ll,2)], [t+st t+st], 'k'); 
       t = t+st;
       xloc = mean([l(ll,1) l(ll,2)]);
       yloc = t;
       if p < 0.01, text(xloc, yloc, '**');
       else, text(xloc, yloc, '*'); end
  end
end
else st = 0; end
ylim([-Inf t+st]);