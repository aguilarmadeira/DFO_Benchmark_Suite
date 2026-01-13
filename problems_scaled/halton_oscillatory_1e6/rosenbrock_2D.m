function varargout = rosenbrock_2D(varargin)
%ROSENBROCK_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   rosenbrock (source instance p=41)
% Dimension: n = 2
% Strategy folder: halton_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 5.12
% Effective contrast: 1996.011964107677
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = rosenbrock_orig(x_orig)

nloc = 2;
lb_orig = [-5.12;-5.12];
ub_orig = [5.12;5.12];
lb_work = [-2562560;-1283.84];
ub_work = [2562560;1283.84];
scale_factors = [500500;250.75];
contrast_ratio = 1996.011964107677;

if nargin == 0
    info.name = mfilename;
    info.problem = 'rosenbrock';
    info.source_p = 41;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;1];
    info.x_global_min_work = [500500;250.75];
    info.global_min_note = 'Rosenbrock (n=2): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = rosenbrock_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = rosenbrock_orig(x);
%
% Purpose:
%
%    Function rosenbrock is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (equal to global minima)
%    Search domain: -1000 <= x <= 1000 (Brachetti et al. (1997))
%                   -2.048 <= x <= 2.048 (Storn and Price (1997))
%                   -5.12 <= x <= 5.12 (Huyer and Neumaier (2008))
%    Cases considered: n = 2 Huyer and Neumaier (2008)
%                            Brachetti et al. (1997)
%                      n = 10 Garcia-Palomares  et al. (2006)
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
n    = size(x,1);
xaux = x(2:n);
f    = sum((x(1:n-1) - 1).^2 + 100.*(x(1:n-1).^2 - xaux).^2,1);
%
% End of rosenbrock.

