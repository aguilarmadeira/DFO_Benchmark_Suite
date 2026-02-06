function varargout = griewank_10D(varargin)
%GRIEWANK_10D  griewank 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (BASELINE HETEROGENEITY):
%
%   x1   ∈ [-400        , 400         ]   (range: 800         )
%   x2   ∈ [-400        , 400         ]   (range: 800         )
%   x3   ∈ [-400        , 400         ]   (range: 800         )
%   x4   ∈ [-400        , 400         ]   (range: 800         )
%   x5   ∈ [-400        , 400         ]   (range: 800         )
%   x6   ∈ [-400        , 400         ]   (range: 800         )
%   x7   ∈ [-400        , 400         ]   (range: 800         )
%   x8   ∈ [-400        , 400         ]   (range: 800         )
%   x9   ∈ [-400        , 400         ]   (range: 800         )
%   x10  ∈ [-400        , 400         ]   (range: 800         )
%
% Effective contrast ratio (max range / min range): 1
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
lb_work = [-400;-400;-400;-400;-400;-400;-400;-400;-400;-400];
ub_work = [400;400;400;400;400;400;400;400;400;400];
scale_factors = [1;1;1;1;1;1;1;1;1;1];
contrast_ratio = 1;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 22;
    info.dimension = nloc;
    info.strategy = 'baseline';
    info.kappa = 1;
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

