function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  schwefel 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-250250000  , 75075000    ]   (range: 325325000   )
%   x2   ∈ [-125375     , 37612.5     ]   (range: 162987.5    )
%   x3   ∈ [-375125000  , 112537500   ]   (range: 487662500   )
%   x4   ∈ [-62937.5    , 18881.25    ]   (range: 81818.75    )
%   x5   ∈ [-312687500  , 93806250    ]   (range: 406493750   )
%   x6   ∈ [-187812.5   , 56343.75    ]   (range: 244156.25   )
%   x7   ∈ [-437562500  , 131268750   ]   (range: 568831250   )
%   x8   ∈ [-31718.75   , 9515.625    ]   (range: 41234.375   )
%   x9   ∈ [-281468750  , 84440625    ]   (range: 365909375   )
%   x10  ∈ [-156593.75  , 46978.125   ]   (range: 203571.875  )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% USAGE:
%   f = schwefel_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = schwefel_10D(n)   % Get bounds for dimension n (must be 10)
%   info = schwefel_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [150;150;150;150;150;150;150;150;150;150];
lb_work = [-250250000;-125375;-375125000;-62937.5;-312687500;-187812.5;-437562500;-31718.75;-281468750;-156593.75];
ub_work = [75075000;37612.5;112537500;18881.25;93806250;56343.75;131268750;9515.625;84440625;46978.125];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5;313.1875];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

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
    info.global_min_known = false;
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

