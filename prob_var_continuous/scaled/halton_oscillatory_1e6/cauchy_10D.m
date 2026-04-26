function varargout = cauchy_10D(varargin)
%CAUCHY_10D  cauchy 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [1001000     , 8808800     ]   (range: 7807800     )
%   x2   ∈ [501.5       , 4413.2      ]   (range: 3911.7      )
%   x3   ∈ [1500500     , 13204400    ]   (range: 11703900    )
%   x4   ∈ [251.75      , 2215.4      ]   (range: 1963.65     )
%   x5   ∈ [1250750     , 11006600    ]   (range: 9755850     )
%   x6   ∈ [751.25      , 6611        ]   (range: 5859.75     )
%   x7   ∈ [1750250     , 15402200    ]   (range: 13651950    )
%   x8   ∈ [126.875     , 1116.5      ]   (range: 989.625     )
%   x9   ∈ [1125875     , 9907700     ]   (range: 8781825     )
%   x10  ∈ [626.375     , 5512.1      ]   (range: 4885.725    )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% USAGE:
%   f = cauchy_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = cauchy_10D(n)   % Get bounds for dimension n (must be 10)
%   info = cauchy_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2;2;2;2;2;2;2;2;2;2];
ub_orig = [17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6];
lb_work = [1001000;501.5000000000002;1500500;251.75;1250750;751.25;1750249.999999999;126.8749999999999;1125875;626.375];
ub_work = [8808800;4413.200000000001;13204400;2215.4;11006600;6611.000000000001;15402200;1116.5;9907700;5512.1];
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
    info.problem = 'cauchy';
    info.source_p = 7;
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

