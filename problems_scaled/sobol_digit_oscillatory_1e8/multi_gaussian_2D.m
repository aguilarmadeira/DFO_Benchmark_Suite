function varargout = multi_gaussian_2D(varargin)
%MULTI_GAUSSIAN_2D  multi_gaussian 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-43992876.5236, 43992876.5236]   (range: 87985753.0472)
%   x2   ∈ [-15945.1744076, 15945.1744076]   (range: 31890.3488152)
%
% Effective contrast ratio (max range / min range): 2759.00880097
%
% Known global minimum (WORK-space):
%   x* = [-298271.7028299645;-108.1082824834812]
%   f* = -1.29695
%
% USAGE:
%   f = multi_gaussian_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = multi_gaussian_2D(n)   % Get bounds for dimension n (must be 2)
%   info = multi_gaussian_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-2;-2];
ub_orig = [2;2];
lb_work = [-43992876.52359302;-15945.17440759307];
ub_work = [43992876.52359302;15945.17440759307];
scale_factors = [21996438.26179651;7972.587203796534];
contrast_ratio = 2759.008800972643;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'multi_gaussian';
    info.source_p = 33;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 2;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1.29695;
    info.x_global_min_orig = [-0.01356;-0.01356];
    info.x_global_min_work = [-298271.7028299645;-108.1082824834812];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Multi-Gaussian (2D): f*=-1.29695. Ref: Ali et al. (2005).';
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

