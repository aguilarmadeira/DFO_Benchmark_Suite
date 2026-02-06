function varargout = wood_4D(varargin)
%WOOD_4D  wood 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [-10         , 10          ]   (range: 20          )
%   x2   ∈ [-10         , 10          ]   (range: 20          )
%   x3   ∈ [-3000       , 3000        ]   (range: 6000        )
%   x4   ∈ [-3000       , 3000        ]   (range: 6000        )
%
% Effective contrast ratio (max range / min range): 300
%
% Known global minimum (WORK-space):
%   x* = [1;1;300.0000000000005;300.0000000000005]
%   f* = 0
%
% USAGE:
%   f = wood_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = wood_4D(n)   % Get bounds for dimension n (must be 4)
%   info = wood_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [10;10;10;10];
lb_work = [-10;-10;-3000;-3000];
ub_work = [10;10;3000;3000];
scale_factors = [1;1;300;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'wood';
    info.source_p = 63;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;1;1;1];
    info.x_global_min_work = [1;1;300.0000000000005;300.0000000000005];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Wood (4D): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 4D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 4 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = wood_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = wood_orig(x);
%
% Purpose:
%
%    Function wood is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [1 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = 100*(x(1,1)^2 - x(2,1))^2 + (1-x(1,1))^2 + 90*(x(3,1)^2-x(4,1))^2 +...
   (1-x(3,1))^2 + 10.1*((1-x(2,1))^2 + (1-x(4,1))^2) + 19.8*(1-x(2,1))*(1-x(4,1));
%
% End of wood.

