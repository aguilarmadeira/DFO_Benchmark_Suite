function varargout = multi_gaussian_2D(varargin)
%MULTI_GAUSSIAN_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   multi_gaussian (source instance p=33)
% Dimension: n = 2
% Strategy folder: halton_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 2
% Effective contrast: 1996.011964107677
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = multi_gaussian_orig(x_orig)

nloc = 2;
lb_orig = [-2;-2];
ub_orig = [2;2];
lb_work = [-1001000;-501.5];
ub_work = [1001000;501.5];
scale_factors = [500500;250.75];
contrast_ratio = 1996.011964107677;

if nargin == 0
    info.name = mfilename;
    info.problem = 'multi_gaussian';
    info.source_p = 33;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 2;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1.29695;
    info.x_global_min_orig = [-0.01356;-0.01356];
    info.x_global_min_work = [-6786.780000000028;-3.400170000000003];
    info.global_min_note = 'Multi-Gaussian (2D): f*=-1.29695. Ref: Ali et al. (2005).';
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
varargout{1} = multi_gaussian_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = multi_gaussian_orig(x);
%
% Purpose:
%
%    Function multi_gaussian is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -1.29695
%    Global minima = [-0.01356 -0.01356]'
%    Local minima (five)
%    Search domain: -2 <= x <= 2 (Montaz et al. (2005))
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
A =[0.5 0 0 0.1
    1.2 1 0 0.5
    1 0 -0.5 0.5
    1 -0.5 0 0.5
    1.2 0 1 0.5];
f = -sum(A(:,1).*exp(-((x(1)-A(:,2)).^2+(x(2)-A(:,3)).^2)./A(:,4).^2),1);
%
% End of multi_gaussian.

