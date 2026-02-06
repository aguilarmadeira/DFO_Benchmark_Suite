function varargout = griewank_10D(varargin)
%GRIEWANK_10D  griewank 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-200200000  , 200200000   ]   (range: 400400000   )
%   x2   ∈ [-100300     , 100300      ]   (range: 200600      )
%   x3   ∈ [-300100000  , 300100000   ]   (range: 600200000   )
%   x4   ∈ [-50350      , 50350       ]   (range: 100700      )
%   x5   ∈ [-250150000  , 250150000   ]   (range: 500300000   )
%   x6   ∈ [-150250     , 150250      ]   (range: 300500      )
%   x7   ∈ [-350050000  , 350050000   ]   (range: 700100000   )
%   x8   ∈ [-25375      , 25375       ]   (range: 50750       )
%   x9   ∈ [-225175000  , 225175000   ]   (range: 450350000   )
%   x10  ∈ [-125275     , 125275      ]   (range: 250550      )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = griewank_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = griewank_10D(n)   % Get bounds for dimension n (must be 10)
%   info = griewank_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-400;-400;-400;-400;-400;-400;-400;-400;-400;-400];
ub_orig = [400;400;400;400;400;400;400;400;400;400];
lb_work = [-200200000;-100300;-300100000;-50350;-250150000;-150250;-350050000;-25375;-225175000;-125275];
ub_work = [200200000;100300;300100000;50350;250150000;150250;350050000;25375;225175000;125275];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5;313.1875];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 22;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 400;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Griewank (n=10): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
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
varargout{1} = griewank_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = griewank_orig(x);
%
% Purpose:
%
%    Function griewank is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -600 <= x <= 600 (Huyer and Neumaier (1999))
%                   -400 <= x <= 400 (Storn and Price (1997))
%    Cases considered: n = 10 Storn and Price (1997)
%                      n = 5 Huyer and Neumaier (1999), ICEO
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
f = 1 + 1/4000*sum(x.^2,1)- prod(cos(x./sqrt([1:n]')),1);
%
% End of griewank.

