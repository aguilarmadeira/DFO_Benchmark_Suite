function varargout = powell_4D(varargin)
%POWELL_4D  powell 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-850178768.918, 850178768.918]   (range: 1700357537.84)
%   x2   ∈ [-42522.1443403, 42522.1443403]   (range: 85044.2886806)
%   x3   ∈ [-57141.8641415, 57141.8641415]   (range: 114283.728283)
%   x4   ∈ [-15325.9096932, 15325.9096932]   (range: 30651.8193863)
%
% Effective contrast ratio (max range / min range): 55473.2988735
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = powell_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = powell_4D(n)   % Get bounds for dimension n (must be 4)
%   info = powell_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [10;10;10;10];
lb_work = [-850178768.917727;-42522.14434031758;-57141.86414149195;-15325.90969317249];
ub_work = [850178768.917727;42522.14434031758;57141.86414149195;15325.90969317249];
scale_factors = [85017876.8917727;4252.214434031758;5714.186414149195;1532.590969317248];
contrast_ratio = 55473.29887350646;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'powell';
    info.source_p = 39;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0];
    info.x_global_min_work = [0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Powell (4D): x*=0, f*=0. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 4D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 4 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = powell_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = powell_orig(x);
%
% Purpose:
%
%    Function powell is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = zeros(4,1)
%    Local minima (equal to global minima)
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
f = (x(1,1) + 10*x(2,1))^2 + 5*(x(3,1) - x(4,1))^2 + (x(2,1) -...
    2*x(3,1))^4 + 10*(x(1,1) - x(4,1))^4;
%
% End of powell.

