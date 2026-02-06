function varargout = fifteenn_local_minima_10D(varargin)
%FIFTEENN_LOCAL_MINIMA_10D  fifteenn_local_minima 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-1492890.625, 1492890.625 ]   (range: 2985781.25  )
%   x2   ∈ [-6487.890625, 6487.890625 ]   (range: 12975.78125 )
%   x3   ∈ [-3990390.625, 3990390.625 ]   (range: 7980781.25  )
%   x4   ∈ [-8985.390625, 8985.390625 ]   (range: 17970.78125 )
%   x5   ∈ [-868515.625 , 868515.625  ]   (range: 1737031.25  )
%   x6   ∈ [-5863.515625, 5863.515625 ]   (range: 11727.03125 )
%   x7   ∈ [-3366015.625, 3366015.625 ]   (range: 6732031.25  )
%   x8   ∈ [-8361.015625, 8361.015625 ]   (range: 16722.03125 )
%   x9   ∈ [-2117265.625, 2117265.625 ]   (range: 4234531.25  )
%   x10  ∈ [-7112.265625, 7112.265625 ]   (range: 14224.53125 )
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   x* = [149289.0625000002;648.7890625000009;399039.0625;898.5390625;86851.56250000012;586.3515625000009;336601.5625000005;836.1015625;211726.5625;711.2265625000009]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = fifteenn_local_minima_10D(n)   % Get bounds for dimension n (must be 10)
%   info = fifteenn_local_minima_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10;10;10];
lb_work = [-1492890.625;-6487.890625;-3990390.625;-8985.390625;-868515.625;-5863.515625;-3366015.625;-8361.015625;-2117265.625;-7112.265625];
ub_work = [1492890.625;6487.890625;3990390.625;8985.390625;868515.625;5863.515625;3366015.625;8361.015625;2117265.625;7112.265625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

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
    info.x_global_min_work = [149289.0625000002;648.7890625000009;399039.0625;898.5390625;86851.56250000012;586.3515625000009;336601.5625000005;836.1015625;211726.5625;711.2265625000009];
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

