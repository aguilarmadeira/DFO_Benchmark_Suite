function varargout = hosaki_2D(varargin)
%HOSAKI_2D  hosaki 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [0           , 5           ]   (range: 5           )
%   x2   ∈ [0           , 1800        ]   (range: 1800        )
%
% Effective contrast ratio (max range / min range): 300
%
% Known global minimum (WORK-space):
%   x* = [4;600]
%   f* = -2.3458
%
% USAGE:
%   f = hosaki_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = hosaki_2D(n)   % Get bounds for dimension n (must be 2)
%   info = hosaki_2D()        % Get complete problem information

nloc = 2;
lb_orig = [0;0];
ub_orig = [5;6];
lb_work = [0;0];
ub_work = [5;1800];
scale_factors = [1;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'hosaki';
    info.source_p = 27;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -2.3458;
    info.x_global_min_orig = [4;2];
    info.x_global_min_work = [4;600];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Hosaki (2D): x*=(4,2), f*=-2.3458. Ref: Ali et al. (2005).';
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
varargout{1} = hosaki_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = hosaki_orig(x);
%
% Purpose:
%
%    Function hosaki is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -2.3458
%    Global minima = [4 2]'
%    Local minima (two)
%    Search domain: 0 <= x(1) <= 5
%                   0 <= x(2) <= 6 (Montaz et al. (2005))
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
f = (1-8*x(1)+7*x(1)^2-7/3*x(1)^3+1/4*x(1)^4)*x(2)^2*exp(-x(2));
%
% End of hosaki.

