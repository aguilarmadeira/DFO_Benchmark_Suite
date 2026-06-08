function varargout = paviani_10D(varargin)
%PAVIANI_10D  paviani 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (PROGRESSIVE HETEROGENEITY):
%
%   x1   ∈ [2.001       , 7.1997      ]   (range: 5.1987      )
%   x2   ∈ [20.01       , 71.997      ]   (range: 51.987      )
%   x3   ∈ [200.1       , 719.97      ]   (range: 519.87      )
%   x4   ∈ [2001        , 7199.7      ]   (range: 5198.7      )
%   x5   ∈ [0.2001      , 0.71997     ]   (range: 0.51987     )
%   x6   ∈ [0.02001     , 0.071997    ]   (range: 0.051987    )
%   x7   ∈ [0.002001    , 0.0071997   ]   (range: 0.0051987   )
%   x8   ∈ [2.001       , 7.1997      ]   (range: 5.1987      )
%   x9   ∈ [20.01       , 71.997      ]   (range: 51.987      )
%   x10  ∈ [200.1       , 719.97      ]   (range: 519.87      )
%
% Effective contrast ratio (max range / min range): 1000000
%
% USAGE:
%   f = paviani_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = paviani_10D(n)   % Get bounds for dimension n (must be 10)
%   info = paviani_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997];
lb_work = [2.000999999999999;20.00999999999999;200.1;2000.999999999999;0.2000999999999999;0.02000999999999999;0.002000999999999999;2.000999999999999;20.00999999999999;200.1];
ub_work = [7.1997;71.997;719.97;7199.7;0.71997;0.07199700000000001;0.0071997;7.1997;71.997;719.97];
scale_factors = [1;10;100;1000;0.1;0.01;0.001;1;10;100];
contrast_ratio = 1000000;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'paviani';
    info.source_p = 36;
    info.dimension = nloc;
    info.strategy = 'progressive';
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
varargout{1} = paviani_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = paviani_orig(x);
%
% Purpose:
%
%    Function paviani is the function described in
%    Montaz et al. (2005).
%
%    dim = 10
%    Minimum global value = -45.778
%    Global minima = 9.351*ones(10,1)
%    Local minima
%    Search domain: 2.001 <= x <= 9.999
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
f = sum((log(x-2)).^2+(log(10-x)).^2,1)-prod(x,1)^0.2;
%
% End of paviani.

