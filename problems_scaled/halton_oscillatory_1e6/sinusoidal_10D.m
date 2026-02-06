function varargout = sinusoidal_10D(varargin)
%SINUSOIDAL_10D  sinusoidal 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 90090000    ]   (range: 90090000    )
%   x2   ∈ [0           , 45135       ]   (range: 45135       )
%   x3   ∈ [0           , 135045000   ]   (range: 135045000   )
%   x4   ∈ [0           , 22657.5     ]   (range: 22657.5     )
%   x5   ∈ [0           , 112567500   ]   (range: 112567500   )
%   x6   ∈ [0           , 67612.5     ]   (range: 67612.5     )
%   x7   ∈ [0           , 157522500   ]   (range: 157522500   )
%   x8   ∈ [0           , 11418.75    ]   (range: 11418.75    )
%   x9   ∈ [0           , 101328750   ]   (range: 101328750   )
%   x10  ∈ [0           , 56373.75    ]   (range: 56373.75    )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% USAGE:
%   f = sinusoidal_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = sinusoidal_10D(n)   % Get bounds for dimension n (must be 10)
%   info = sinusoidal_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [180;180;180;180;180;180;180;180;180;180];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [90090000;45135;135045000;22657.5;112567500;67612.5;157522500;11418.75;101328750;56373.75];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5;313.1875];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'sinusoidal';
    info.source_p = 53;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
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
varargout{1} = sinusoidal_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = sinusoidal_orig(x);
%
% Purpose:
%
%    Function sinusoidal is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -(A+1)
%    Global minima = (90+z)*ones(n,1)
%    Local minima (several)
%    Search domain: 0 <= x <= 180 (Montaz et al. (2005))
%    Cases considered: n = 10, 20
%                      A = 2.5;B = 5;z = 30;Montaz et al. (2005)
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
A = 2.5;
B = 5;
z = 30;
f = -(A*prod(sind(x-z),1)+prod(sind(B.*(x-z)),1));
%
% End of sinusoidal.

