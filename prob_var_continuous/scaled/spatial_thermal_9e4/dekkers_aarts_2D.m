function varargout = dekkers_aarts_2D(varargin)
%DEKKERS_AARTS_2D  dekkers_aarts 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [-20         , 6           ]   (range: 26          )
%   x2   ∈ [-6000       , 1800        ]   (range: 7800        )
%
% Effective contrast ratio (max range / min range): 300
%
% USAGE:
%   f = dekkers_aarts_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = dekkers_aarts_2D(n)   % Get bounds for dimension n (must be 2)
%   info = dekkers_aarts_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-20;-20];
ub_orig = [6;6];
lb_work = [-20;-6000];
ub_work = [6;1800];
scale_factors = [1;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'dekkers_aarts';
    info.source_p = 10;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 20;  % Original bound(p) from test setup
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
varargout{1} = dekkers_aarts_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = dekkers_aarts_orig(x);
%
% Purpose:
%
%    Function dekkers_aarts is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -24771
%    Global minima = [0 15]', [0 -15]'
%    Local minima (three)
%    Search domain: -20 <= x <= 20 (Montaz et al. (2005))
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
f = 10^5*x(1)^2+x(2)^2-(x(1)^2+x(2)^2)^2+1/10^5*(x(1)^2+x(2)^2)^4;
%
% End of dekkers_aarts.

