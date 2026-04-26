function varargout = cauchy_10D(varargin)
%CAUCHY_10D  cauchy 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (PROGRESSIVE HETEROGENEITY):
%
%   x1   ∈ [2           , 17.6        ]   (range: 15.6        )
%   x2   ∈ [20          , 176         ]   (range: 156         )
%   x3   ∈ [200         , 1760        ]   (range: 1560        )
%   x4   ∈ [2000        , 17600       ]   (range: 15600       )
%   x5   ∈ [0.2         , 1.76        ]   (range: 1.56        )
%   x6   ∈ [0.02        , 0.176       ]   (range: 0.156       )
%   x7   ∈ [0.002       , 0.0176      ]   (range: 0.0156      )
%   x8   ∈ [2           , 17.6        ]   (range: 15.6        )
%   x9   ∈ [20          , 176         ]   (range: 156         )
%   x10  ∈ [200         , 1760        ]   (range: 1560        )
%
% Effective contrast ratio (max range / min range): 1000000
%
% USAGE:
%   f = cauchy_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = cauchy_10D(n)   % Get bounds for dimension n (must be 10)
%   info = cauchy_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2;2;2;2;2;2;2;2;2;2];
ub_orig = [17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6];
lb_work = [2;20;200;1999.999999999999;0.2;0.01999999999999999;0.002000000000000001;2;20;200];
ub_work = [17.6;176;1760;17600;1.76;0.176;0.0176;17.6;176;1760];
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
    info.problem = 'cauchy';
    info.source_p = 7;
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
varargout{1} = cauchy_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = cauchy_orig(x);
%
% Purpose:
%
%    Function cauchy is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value
%    Global minima
%    Local minima
%    Search domain: y(1) <= x <= y(dim) (Brachetti et al. (1997))
%    Cases considered: n = 4 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
%                      n = 25 Brachetti et al. (1997)
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
n = size(x,1);
if n == 4
    y = [3 7 12 17]';
else
    if n == 10
        y = [2 5 7 8 11 15 17 21 23 26]';
    else
        if n == 25
            y = [4.1 7.7 17.5 31.4 32.7 92.4 115.3 118.3 119 129.6 198.6];
            y = [y, 200.7 242.5 255 274.7 274.7 303.8 334.1 430 489.1 703.4 ];
            y = [y, 978 1656 1697.8 2745.6]';
        end
    end
end
f = -sum(log(pi)+log(1+(y-x).^2),1);
%
% End of cauchy.

