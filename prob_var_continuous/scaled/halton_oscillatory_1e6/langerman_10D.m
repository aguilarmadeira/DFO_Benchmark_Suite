function varargout = langerman_10D(varargin)
%LANGERMAN_10D  langerman 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 3253250     ]   (range: 3253250     )
%   x2   ∈ [0           , 1629.875    ]   (range: 1629.875    )
%   x3   ∈ [0           , 4876625     ]   (range: 4876625     )
%   x4   ∈ [0           , 818.1875    ]   (range: 818.1875    )
%   x5   ∈ [0           , 4064937.5   ]   (range: 4064937.5   )
%   x6   ∈ [0           , 2441.5625   ]   (range: 2441.5625   )
%   x7   ∈ [0           , 5688312.5   ]   (range: 5688312.5   )
%   x8   ∈ [0           , 412.34375   ]   (range: 412.34375   )
%   x9   ∈ [0           , 3659093.75  ]   (range: 3659093.75  )
%   x10  ∈ [0           , 2035.71875  ]   (range: 2035.71875  )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% USAGE:
%   f = langerman_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = langerman_10D(n)   % Get bounds for dimension n (must be 10)
%   info = langerman_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [6.5;6.5;6.5;6.5;6.5;6.5;6.5;6.5;6.5;6.5];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [3253250;1629.875;4876625;818.1875;4064937.5;2441.5625;5688312.5;412.34375;3659093.75;2035.71875];
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
    info.problem = 'langerman';
    info.source_p = 29;
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
varargout{1} = langerman_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = langerman_orig(x);
%
% Purpose:
%
%    Function langerman is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -0.965, for n = 5, 10
%    Global minima: [8.074 8.777 3.467 1.867 6.708]' for n = 5
%                   [8.074 8.777 3.467 1.867 6.708 6.349 4.534...
%                    0.276 7.633 1.567]' for n = 10
%    Local minima 
%    Search domain: 0 <= x <= 10 (Montaz et al. (2005))
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
dim = size(x,1);
A   = [0.806 9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020
0.517 9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374
0.100 8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982
0.908 2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426
0.965 8.074 8.777 3.467 1.867 6.708 6.349 4.534 0.276 7.633 1.567];
d = sum((repmat(x',5,1)-A(:,2:11)).^2,2);
f = -sum(A(:,1).*cos(pi.*d).*exp(-d./pi),1);
%
% End of langerman.

