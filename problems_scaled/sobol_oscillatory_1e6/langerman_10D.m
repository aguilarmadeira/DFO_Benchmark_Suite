function varargout = langerman_10D(varargin)
%LANGERMAN_10D  langerman 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 1492890.625 ]   (range: 1492890.625 )
%   x2   ∈ [0           , 6487.890625 ]   (range: 6487.890625 )
%   x3   ∈ [0           , 3990390.625 ]   (range: 3990390.625 )
%   x4   ∈ [0           , 8985.390625 ]   (range: 8985.390625 )
%   x5   ∈ [0           , 868515.625  ]   (range: 868515.625  )
%   x6   ∈ [0           , 5863.515625 ]   (range: 5863.515625 )
%   x7   ∈ [0           , 3366015.625 ]   (range: 3366015.625 )
%   x8   ∈ [0           , 8361.015625 ]   (range: 8361.015625 )
%   x9   ∈ [0           , 2117265.625 ]   (range: 2117265.625 )
%   x10  ∈ [0           , 7112.265625 ]   (range: 7112.265625 )
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   x* = [1205359.890625;5694.4216015625;1383468.4296875;1677.5724296875;582600.28125;3722.7460703125;1526151.484375;230.76403125;1616108.8515625;1114.4920234375]
%   f* = -0.965
%
% USAGE:
%   f = langerman_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = langerman_10D(n)   % Get bounds for dimension n (must be 10)
%   info = langerman_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [10;10;10;10;10;10;10;10;10;10];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [1492890.625;6487.890625;3990390.625;8985.390625;868515.625;5863.515625;3366015.625;8361.015625;2117265.625;7112.265625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'langerman';
    info.source_p = 29;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -0.965;
    info.x_global_min_orig = [8.074;8.776999999999999;3.467;1.867;6.708;6.349;4.534;0.276;7.633;1.567];
    info.x_global_min_work = [1205359.890625;5694.4216015625;1383468.4296875;1677.5724296875;582600.28125;3722.7460703125;1526151.484375;230.76403125;1616108.8515625;1114.4920234375];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Langerman (10D): f*=-0.965. Ref: Ali et al. (2005).';
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
varargout{1} = langerman_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = langerman_orig(x);
%
% Purpose:
%
%    Function langerman is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -0.965, for n = 5, 10
%    Global minima: [8.074 8.777 3.467 1.867 6.708]' for n = 5
%                   [8.074 8.777 3.467 1.867 6.708 6.349 4.534...
%                    0.276 7.633 1.567]' for n = 10
%    Local minima 
%    Search domain: 0 <= x <= 10 (Montaz et al. (2005))
%    Cases considered: n = 10 Montaz et al. (2005)
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
A   = [0.806 9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020
0.517 9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374
0.100 8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982
0.908 2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426
0.965 8.074 8.777 3.467 1.867 6.708 6.349 4.534 0.276 7.633 1.567];
d = sum((repmat(x',5,1)-A(:,2:11)).^2,2);
f = -sum(A(:,1).*cos(pi.*d).*exp(-d./pi),1);
%
% End of langerman.

