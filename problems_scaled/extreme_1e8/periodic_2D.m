function varargout = periodic_2D(varargin)
%PERIODIC_2D  periodic 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (EXTREME HETEROGENEITY):
%
%   x1   ∈ [-10         , 10          ]   (range: 20          )
%   x2   ∈ [-1000000000 , 1000000000  ]   (range: 2000000000  )
%
% Effective contrast ratio (max range / min range): 100000000
%
% Known global minimum (WORK-space):
%   x* = [0;0]
%   f* = 0.9
%
% USAGE:
%   f = periodic_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = periodic_2D(n)   % Get bounds for dimension n (must be 2)
%   info = periodic_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-10;-10];
ub_orig = [10;10];
lb_work = [-10;-1000000000];
ub_work = [10;1000000000];
scale_factors = [1;100000000];
contrast_ratio = 100000000;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'periodic';
    info.source_p = 37;
    info.dimension = nloc;
    info.strategy = 'extreme';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0.9;
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Periodic (2D): x*=0, f*=0.9. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 2D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 2 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = periodic_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = periodic_orig(x);
%
% Purpose:
%
%    Function periodic is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0.9
%    Global minima = [0 0]'
%    Local minima (49 with value 1)
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
f = 1+sin(x(1))^2+sin(x(2))^2-0.1*exp(-x(1)^2-x(2)^2);
%
% End of periodic.

