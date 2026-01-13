function varargout = powell_4D(varargin)
%POWELL_4D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   powell (source instance p=39)
% Dimension: n = 4
% Strategy folder: progressive (kappa = 1000000)
% Original bound tag: bound(p) = 10
% Effective contrast: 1000
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = powell_orig(x_orig)

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [10;10;10;10];
lb_work = [-10;-100;-1000;-10000];
ub_work = [10;100;1000;10000];
scale_factors = [1;10;100;1000];
contrast_ratio = 1000;

if nargin == 0
    info.name = mfilename;
    info.problem = 'powell';
    info.source_p = 39;
    info.dimension = nloc;
    info.strategy = 'progressive';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0];
    info.x_global_min_work = [0;0;0;0];
    info.global_min_note = 'Powell (4D): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = powell_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = powell_orig(x);
%
% Purpose:
%
%    Function powell is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = zeros(4,1)
%    Local minima (equal to global minima)
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
f = (x(1,1) + 10*x(2,1))^2 + 5*(x(3,1) - x(4,1))^2 + (x(2,1) -...
    2*x(3,1))^4 + 10*(x(1,1) - x(4,1))^4;
%
% End of powell.

