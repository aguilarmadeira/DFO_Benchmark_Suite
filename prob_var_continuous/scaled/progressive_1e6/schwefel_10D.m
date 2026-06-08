function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  schwefel 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (PROGRESSIVE HETEROGENEITY):
%
%   x1   ∈ [-500        , 150         ]   (range: 650         )
%   x2   ∈ [-5000       , 1500        ]   (range: 6500        )
%   x3   ∈ [-50000      , 15000       ]   (range: 65000       )
%   x4   ∈ [-500000     , 150000      ]   (range: 650000      )
%   x5   ∈ [-50         , 15          ]   (range: 65          )
%   x6   ∈ [-5          , 1.5         ]   (range: 6.5         )
%   x7   ∈ [-0.5        , 0.15        ]   (range: 0.65        )
%   x8   ∈ [-500        , 150         ]   (range: 650         )
%   x9   ∈ [-5000       , 1500        ]   (range: 6500        )
%   x10  ∈ [-50000      , 15000       ]   (range: 65000       )
%
% Effective contrast ratio (max range / min range): 1000000
%
% USAGE:
%   f = schwefel_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = schwefel_10D(n)   % Get bounds for dimension n (must be 10)
%   info = schwefel_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [150;150;150;150;150;150;150;150;150;150];
lb_work = [-500;-5000;-50000;-500000;-50;-5;-0.5;-500;-5000;-50000];
ub_work = [150;1500;15000;150000;15;1.5;0.15;150;1500;15000];
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
    info.problem = 'schwefel';
    info.source_p = 46;
    info.dimension = nloc;
    info.strategy = 'progressive';
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

