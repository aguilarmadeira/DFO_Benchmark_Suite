function varargout = aluffi_pentini_2D(varargin)
%ALUFFI_PENTINI_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   aluffi_pentini (source instance p=2)
% Dimension: n = 2
% Strategy folder: progressive (kappa = 1000000)
% Original bound tag: bound(p) = 10
% Effective contrast: 10
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = aluffi_pentini_orig(x_orig)

nloc = 2;
lb_orig = [-10;-10];
ub_orig = [10;10];
lb_work = [-10;-100];
ub_work = [10;100];
scale_factors = [1;10];
contrast_ratio = 10;

if nargin == 0
    info.name = mfilename;
    info.problem = 'aluffi_pentini';
    info.source_p = 2;
    info.dimension = nloc;
    info.strategy = 'progressive';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -0.3523;
    info.x_global_min_orig = [-1.0465;0];
    info.x_global_min_work = [-1.0465;0];
    info.global_min_note = 'Aluffi-Pentini (2D): x*=(-1.0465,0), f*=-0.3523. Ref: Ali et al. (2005).';
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
varargout{1} = aluffi_pentini_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = aluffi_pentini_orig(x);
%
% Purpose:
%
%    Function aluffi_pentini is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -0.3523
%    Global minima = [-1.0465 0]'
%    Local minima (two)
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
f = 0.25*x(1)^4-0.5*x(1)^2+0.1*x(1)+0.5*x(2)^2;
%
% End of aluffi_pentini.

