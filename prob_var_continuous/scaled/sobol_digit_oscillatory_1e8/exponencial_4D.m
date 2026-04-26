function varargout = exponencial_4D(varargin)
%EXPONENCIAL_4D  exponencial 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-14552626.7315, 4365788.01944]   (range: 18918414.7509)
%   x2   ∈ [-7120.79244681, 2136.23773404]   (range: 9257.03018085)
%   x3   ∈ [-3700.79769291, 1110.23930787]   (range: 4811.03700078)
%   x4   ∈ [-8038.83172702, 2411.6495181]   (range: 10450.4812451)
%
% Effective contrast ratio (max range / min range): 3932.29458593
%
% Known global minimum (WORK-space):
%   x* = [-1.862645149230957e-09;0;-4.547473508864641e-13;0]
%   f* = -1
%
% USAGE:
%   f = exponencial_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = exponencial_4D(n)   % Get bounds for dimension n (must be 4)
%   info = exponencial_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-1;-1;-1;-1];
ub_orig = [0.3;0.3;0.3;0.3];
lb_work = [-14552626.73146438;-7120.792446808101;-3700.797692908849;-8038.831727015808];
ub_work = [4365788.019439314;2136.237734042431;1110.239307872655;2411.649518104743];
scale_factors = [14552626.73146438;7120.792446808101;3700.797692908849;8038.831727015808];
contrast_ratio = 3932.294585934506;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'exponencial';
    info.source_p = 14;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 1;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1;
    info.x_global_min_orig = [0;0;0;0];
    info.x_global_min_work = [-1.862645149230957e-09;0;-4.547473508864641e-13;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Exponencial (n=4): x*=0, f*=-1. Ref: Brachetti et al. (1997).';
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
varargout{1} = exponencial_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = exponencial_orig(x);
%
% Purpose:
%
%    Function exponencial is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = -1
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -1 <= x <= 1 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
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
f = -exp(-0.5*sum(x.^2,1)); 
%
% End of exponencial.

