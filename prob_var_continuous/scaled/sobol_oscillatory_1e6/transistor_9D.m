function varargout = transistor_9D(varargin)
%TRANSISTOR_9D  transistor 9D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-1492890.625, 447867.1875 ]   (range: 1940757.8125)
%   x2   ∈ [-6487.890625, 1946.3671875]   (range: 8434.2578125)
%   x3   ∈ [-3990390.625, 1197117.1875]   (range: 5187507.8125)
%   x4   ∈ [-8985.390625, 2695.6171875]   (range: 11681.0078125)
%   x5   ∈ [-868515.625 , 260554.6875 ]   (range: 1129070.3125)
%   x6   ∈ [-5863.515625, 1759.0546875]   (range: 7622.5703125)
%   x7   ∈ [-3366015.625, 1009804.6875]   (range: 4375820.3125)
%   x8   ∈ [-8361.015625, 2508.3046875]   (range: 10869.3203125)
%   x9   ∈ [-2117265.625, 635179.6875 ]   (range: 2752445.3125)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% USAGE:
%   f = transistor_9D(x)          % Evaluate function at point x (9D vector)
%   [lb, ub] = transistor_9D(n)   % Get bounds for dimension n (must be 9)
%   info = transistor_9D()        % Get complete problem information

nloc = 9;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3;3;3;3];
lb_work = [-1492890.625;-6487.890625;-3990390.625;-8985.390625;-868515.625;-5863.515625;-3366015.625;-8361.015625;-2117265.625];
ub_work = [447867.1875;1946.3671875;1197117.1875;2695.6171875;260554.6875;1759.0546875;1009804.6875;2508.3046875;635179.6875];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'transistor';
    info.source_p = 62;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
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
    if arg1 ~= nloc, error('This instance is 9D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 9 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = transistor_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = transistor_orig(x);
%
% Purpose:
%
%    Function transistor is the function described in
%    Montaz et al. (2005).
%
%    dim = 9
%    Minimum global value = 0
%    Global minima near [0.9 0.45 1 2 8 8 5 1 2]'
%    Local minima
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
g = [0.485 0.752 0.869 0.982
    0.369 1.254 0.703 1.455
    5.2095 10.0677 22.9274 20.2153
    23.3037 101.779 111.461 191.267
    28.5132 111.8467 134.3884 211.4823]';
alfa = (1-x(1)*x(2))*x(3).*(exp(x(5).*(g(:,1)-g(:,3).*x(7)*0.001-...
       g(:,5).*x(8)*0.001))-1)-g(:,5)+g(:,4).*x(2);
beta = (1-x(1)*x(2))*x(4).*(exp(x(6).*(g(:,1)-g(:,2)-g(:,3).*x(7)*0.001+...
       g(:,4).*x(9)*0.001))-1)-g(:,5).*x(1)+g(:,4);
f = (x(1)*x(3)-x(2)*x(4))^2 + sum(alfa.^2+beta.^2);
%
% End of transistor.

