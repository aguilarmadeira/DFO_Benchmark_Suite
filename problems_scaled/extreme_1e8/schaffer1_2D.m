function varargout = schaffer1_2D(varargin)
%SCHAFFER1_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   schaffer1 (source instance p=44)
% Dimension: n = 2
% Strategy folder: extreme (kappa = 100000000)
% Original bound tag: bound(p) = 100
% Effective contrast: 100000000
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = schaffer1_orig(x_orig)

nloc = 2;
lb_orig = [-100;-100];
ub_orig = [100;100];
lb_work = [-100;-10000000000];
ub_work = [100;10000000000];
scale_factors = [1;100000000];
contrast_ratio = 100000000;

if nargin == 0
    info.name = mfilename;
    info.problem = 'schaffer1';
    info.source_p = 44;
    info.dimension = nloc;
    info.strategy = 'extreme';
    info.kappa = 100000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;0];
    info.global_min_note = 'Schaffer1 (2D): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = schaffer1_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schaffer1_orig(x);
%
% Purpose:
%
%    Function schaffer1 is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
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
f = 0.5 + (sin(norm(x,2))^2-0.5)/((1+0.001*(x(1)^2+x(2)^2))^2);
%
% End of schaffer1.

