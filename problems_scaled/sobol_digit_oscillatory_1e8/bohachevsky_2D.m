function varargout = bohachevsky_2D(varargin)
%BOHACHEVSKY_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   bohachevsky (source instance p=4)
% Dimension: n = 2
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 50
% Effective contrast: 2.225017544428539
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = bohachevsky_orig(x_orig)

nloc = 2;
lb_orig = [-50;-50];
ub_orig = [50;50];
lb_work = [-214561.3712491645;-477402.8153860362];
ub_work = [214561.3712491645;477402.8153860362];
scale_factors = [4291.22742498329;9548.056307720724];
contrast_ratio = 2.225017544428539;

if nargin == 0
    info.name = mfilename;
    info.problem = 'bohachevsky';
    info.source_p = 4;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 50;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;0];
    info.global_min_note = 'Bohachevsky (2D): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = bohachevsky_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = bohachevsky_orig(x);
%
% Purpose:
%
%    Function bohachevsky is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -50 <= x <= 50 (Montaz et al. (2005))
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
f = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2))+0.3;
%
% End of bohachevsky.

