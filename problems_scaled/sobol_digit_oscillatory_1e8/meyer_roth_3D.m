function varargout = meyer_roth_3D(varargin)
%MEYER_ROTH_3D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   meyer_roth (source instance p=31)
% Dimension: n = 3
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 10
% Effective contrast: 214178.3388140479
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = meyer_roth_orig(x_orig)

nloc = 3;
lb_orig = [-10;-10;-10];
ub_orig = [10;10;10];
lb_work = [-876263067.9601548;-4091.2777305689;-60349.83185248793];
ub_work = [876263067.9601548;4091.2777305689;60349.83185248793];
scale_factors = [87626306.79601547;409.12777305689;6034.983185248793];
contrast_ratio = 214178.3388140479;

if nargin == 0
    info.name = mfilename;
    info.problem = 'meyer_roth';
    info.source_p = 31;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
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
    if arg1 ~= nloc, error('This instance is 3D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 3 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = meyer_roth_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = meyer_roth_orig(x);
%
% Purpose:
%
%    Function meyer_roth is the function described in
%    Brachetti et al. (1997).
%
%    dim = 3
%    Minimum global value = 0.4 * 10^-4
%    Global minima = [3.13 15.16 0.78]'
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
T = [1 2 1 2 0.1]';
V = [1 1 2 2 0]';
Y = [0.126 0.219 0.076 0.126 0.186]';
f = sum(((x(1)*x(3).*T)./(1+x(1).*T+x(2).*V)-Y).^2,1);
%
% End of meyer_roth.

