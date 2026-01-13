function varargout = shubert_2D(varargin)
%SHUBERT_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   shubert (source instance p=52)
% Dimension: n = 2
% Strategy folder: baseline (kappa = 1)
% Original bound tag: bound(p) = 10
% Effective contrast: 1
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = shubert_orig(x_orig)

nloc = 2;
lb_orig = [-10;-10];
ub_orig = [10;10];
lb_work = [-10;-10];
ub_work = [10;10];
scale_factors = [1;1];
contrast_ratio = 1;

if nargin == 0
    info.name = mfilename;
    info.problem = 'shubert';
    info.source_p = 52;
    info.dimension = nloc;
    info.strategy = 'baseline';
    info.kappa = 1;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -186.7309;
    info.global_min_note = 'Shubert (2D): f*=-186.7309, 18 global minima. Ref: Ali et al. (2005).';
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
varargout{1} = shubert_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = shubert_orig(x);
%
% Purpose:
%
%    Function shubert is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -186.7309
%    Global minima (18)
%    Local minima (760)
%    Search domain: -10 <= x <= 10 (Huyer and Neumaier (1999))
%                                  (Huyer and Neumaier (2008))
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
f = sum([1:5].*cos([2:6].*x(1)+[1:5]),2)*sum([1:5].*cos([2:6].*x(2)+[1:5]),2);
%
% End of shubert.

