function varargout = poissonian_2D(varargin)
%POISSONIAN_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   poissonian (source instance p=38)
% Dimension: n = 2
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 1.705911702308682
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = poissonian_orig(x_orig)

nloc = 2;
lb_orig = [1;1];
ub_orig = [21;8];
lb_work = [76779010.2413801;45007610.96689337];
ub_work = [1612359215.06898;360060887.7351471];
scale_factors = [76779010.24138001;45007610.96689338];
contrast_ratio = 1.705911702308682;

if nargin == 0
    info.name = mfilename;
    info.problem = 'poissonian';
    info.source_p = 38;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -95.28;
    info.global_min_note = 'Poissonian (2D): f*=-95.28, x* not documented. Ref: Brachetti et al. (1997).';
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
varargout{1} = poissonian_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = poissonian_orig(x);
%
% Purpose:
%
%    Function poissonian is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = -95.28
%    Global minima
%    Local minima
%    Search domain: 1 <= x(1) <= 21 
%                   1 <= x(2) <= 8 (Brachetti et al. (1997))
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
n      = [5 2 4 2 7 2 4 5 4 4 15 10 8 15 5 6 3 4 5 2 6]';
i      = [1:21]';
lambda = 2+5.*exp(-0.5.*(((i-x(1))./x(2)).^2))+3;
f      = sum(-lambda + n.*log(lambda),1);
%
% End of poissonian.

