function varargout = mccormick_2D(varargin)
%MCCORMICK_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   mccormick (source instance p=30)
% Dimension: n = 2
% Strategy folder: halton_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 1996.011964107677
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = mccormick_orig(x_orig)

nloc = 2;
lb_orig = [-1.5;-3];
ub_orig = [4;3];
lb_work = [-750750;-752.25];
ub_work = [2002000;752.25];
scale_factors = [500500;250.75];
contrast_ratio = 1996.011964107677;

if nargin == 0
    info.name = mfilename;
    info.problem = 'mccormick';
    info.source_p = 30;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1.9133;
    info.x_global_min_orig = [-0.547;-1.547];
    info.x_global_min_work = [-273773.5;-387.91025];
    info.global_min_note = 'McCormick (2D): f*=-1.9133. Ref: Ali et al. (2005).';
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
varargout{1} = mccormick_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = mccormick_orig(x);
%
% Purpose:
%
%    Function mccormick is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -1.9133
%    Global minima = [-0.547 -1.547]'
%    Local minima (two)
%    Search domain: -1.5 <= x(1) <= 4
%                     -3 <= x(2) <= 3 (Montaz et al. (2005))
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
f = sin(x(1)+x(2))+(x(1)-x(2))^2-3/2*x(1)+5/2*x(2)+1;
%
% End of mccormick.

