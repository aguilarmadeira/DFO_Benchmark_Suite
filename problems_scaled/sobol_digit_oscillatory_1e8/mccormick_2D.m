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
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 26382.20497486475
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = mccormick_orig(x_orig)

nloc = 2;
lb_orig = [-1.5;-3];
ub_orig = [4;3];
lb_work = [-4745.523290163146;-250394736.3081574];
ub_work = [12654.72877376839;250394736.3081574];
scale_factors = [3163.682193442098;83464912.10271913];
contrast_ratio = 26382.20497486475;

if nargin == 0
    info.name = mfilename;
    info.problem = 'mccormick';
    info.source_p = 30;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1.9133;
    info.x_global_min_orig = [-0.547;-1.547];
    info.x_global_min_work = [-1730.534159812827;-129120219.0229065];
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

