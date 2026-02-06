function varargout = griewank_10D(varargin)
%GRIEWANK_10D  griewank 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-59715625   , 59715625    ]   (range: 119431250   )
%   x2   ∈ [-259515.625 , 259515.625  ]   (range: 519031.25   )
%   x3   ∈ [-159615625  , 159615625   ]   (range: 319231250   )
%   x4   ∈ [-359415.625 , 359415.625  ]   (range: 718831.25   )
%   x5   ∈ [-34740625   , 34740625    ]   (range: 69481250    )
%   x6   ∈ [-234540.625 , 234540.625  ]   (range: 469081.25   )
%   x7   ∈ [-134640625  , 134640625   ]   (range: 269281250   )
%   x8   ∈ [-334440.625 , 334440.625  ]   (range: 668881.25   )
%   x9   ∈ [-84690625   , 84690625    ]   (range: 169381250   )
%   x10  ∈ [-284490.625 , 284490.625  ]   (range: 568981.25   )
%
% Effective contrast ratio (max range / min range): 680.545747672
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
lb_work = [-59715625;-259515.625;-159615625;-359415.625;-34740625;-234540.625;-134640625;-334440.625;-84690625;-284490.625];
ub_work = [59715625;259515.625;159615625;359415.625;34740625;234540.625;134640625;334440.625;84690625;284490.625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 22;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
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

