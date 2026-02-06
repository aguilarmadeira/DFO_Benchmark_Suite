function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  schwefel 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-250250000  , 250250000   ]   (range: 500500000   )
%   x2   ∈ [-125375     , 125375      ]   (range: 250750      )
%   x3   ∈ [-375125000  , 375125000   ]   (range: 750250000   )
%   x4   ∈ [-62937.5    , 62937.5     ]   (range: 125875      )
%   x5   ∈ [-312687500  , 312687500   ]   (range: 625375000   )
%   x6   ∈ [-187812.5   , 187812.5    ]   (range: 375625      )
%   x7   ∈ [-437562500  , 437562500   ]   (range: 875125000   )
%   x8   ∈ [-31718.75   , 31718.75    ]   (range: 63437.5     )
%   x9   ∈ [-281468750  , 281468750   ]   (range: 562937500   )
%   x10  ∈ [-156593.75  , 156593.75   ]   (range: 313187.5    )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [210694834.35;105557.901525;315831767.175;52989.4351125;263263300.7624999;158126.3679375;368400233.5875;26705.20190625;236979067.55625;131842.13473125]
%   f* = -4189.829
%
% USAGE:
%   f = schwefel_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = schwefel_10D(n)   % Get bounds for dimension n (must be 10)
%   info = schwefel_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [500;500;500;500;500;500;500;500;500;500];
lb_work = [-250250000;-125375;-375125000;-62937.5;-312687500;-187812.5;-437562500;-31718.75;-281468750;-156593.75];
ub_work = [250250000;125375;375125000;62937.5;312687500;187812.5;437562500;31718.75;281468750;156593.75];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5;313.1875];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'schwefel';
    info.source_p = 46;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 500;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -4189.829;
    info.x_global_min_orig = [420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687];
    info.x_global_min_work = [210694834.35;105557.901525;315831767.175;52989.4351125;263263300.7624999;158126.3679375;368400233.5875;26705.20190625;236979067.55625;131842.13473125];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Schwefel (n=10): x*=420.9687, f*=-418.9829*n. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 10D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 10 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = schwefel_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schwefel_orig(x);
%
% Purpose:
%
%    Function schwefel is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -418.9829*n
%    Global minima = 420.97*ones(n,1)
%    Local minima 
%    Search domain: -500 <= x <= 500 (Montaz et al. (2005))
%    Cases considered: n = 10 Montaz et al. (2005)
%
% Input:  
%
%         x (point given by the optimizer).
%
% Output: 
%
%         f (function value at x).
%
% Written by A. L. Custodio and J. F. A. Madeira.
%
% Version June 2012.
%
%
f = -sum(x.*sin(sqrt(abs(x))),1);
%
% End of schwefel.

