function varargout = fifteenn_local_minima_10D(varargin)
%FIFTEENN_LOCAL_MINIMA_10D  fifteenn_local_minima 10D test problem (heterogeneous WORK-space wrapper).
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
%   x10  ∈ [-7112.265625, 2133.6796875]   (range: 9245.9453125)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   x* = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = fifteenn_local_minima_10D(n)   % Get bounds for dimension n (must be 10)
%   info = fifteenn_local_minima_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3;3;3;3;3];
lb_work = [-1492890.625;-6487.890625;-3990390.625;-8985.390625;-868515.625;-5863.515625;-3366015.625;-8361.015625;-2117265.625;-7112.265625];
ub_work = [447867.1875;1946.3671875;1197117.1875;2695.6171875;260554.6875;1759.0546875;1009804.6875;2508.3046875;635179.6875;2133.6796875];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 19;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;1;1;1;1;1;1;1;1;1];
    info.x_global_min_work = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=10): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 10D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 10 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = fifteenn_local_minima_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = fifteenn_local_minima_orig(x);
%
% Purpose:
%
%    Function fifteenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (15^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
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
dim = size(x,1);
f   = 0.1*sum(((x(1:dim-1)-1).^2).*(1+10*sin(3*pi.*x(2:dim)).^2),1);
f   = f + 0.1*(sin(3*pi*x(1))^2+(x(dim)-1)^2*(1+sin(2*pi*x(dim))^2));
%
% End of fifteenn_local_minima.

