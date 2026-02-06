function varargout = branin_hoo_2D(varargin)
%BRANIN_HOO_2D  branin_hoo 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [-5          , 10          ]   (range: 15          )
%   x2   ∈ [0           , 4500        ]   (range: 4500        )
%
% Effective contrast ratio (max range / min range): 300
%
% Known global minimum (WORK-space):
%   multiple minimizers; see info.x_global_min_work
%   f* = 0.3978873577297384
%
% USAGE:
%   f = branin_hoo_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = branin_hoo_2D(n)   % Get bounds for dimension n (must be 2)
%   info = branin_hoo_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-5;0];
ub_orig = [10;15];
lb_work = [-5;0];
ub_work = [10;4500];
scale_factors = [1;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'branin_hoo';
    info.source_p = 5;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0.3978873577297384;
    info.x_global_min_orig = [-3.141592653589793 3.141592653589793 9.424777960769379;12.275 2.275 2.475];
    info.x_global_min_work = [-3.141592653589793 3.141592653589793 9.424777960769379;3682.5 682.5 742.5];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Branin-Hoo (2D): 3 global minima, f*=5/(4*pi). Ref: Garcia-Palomares (2012).';
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
varargout{1} = branin_hoo_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = branin_hoo_orig(x);
%
% Purpose:
%
%    Function brannin_hoo is the function described in
%    Garcia-Palomares(2012).
%
%    dim = 2
%    Minimum global value =  5/(4*pi)
%    Global minima = ([-pi 12.275]', [pi 2.275]', [3*pi 2.475]')
%    Local minima (equal to global minima)
%    Search domain: -5 <= x <= 20 (Garcia-Palomares(2012))
%
%                    -5 <= x(1) <= 10 
%                     0 <= x(2) <= 15 (Huyer and Neumaier (1999))
%                                     (Huyer and Neumaier (2008))
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
f = (x(2)-1.275*(x(1)/pi)^2+5*x(1)/pi-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
%
% End of branin_hoo.

