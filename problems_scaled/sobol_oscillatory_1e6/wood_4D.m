function varargout = wood_4D(varargin)
%WOOD_4D  wood 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-1492890.625, 1492890.625 ]   (range: 2985781.25  )
%   x2   ∈ [-6487.890625, 6487.890625 ]   (range: 12975.78125 )
%   x3   ∈ [-3990390.625, 3990390.625 ]   (range: 7980781.25  )
%   x4   ∈ [-8985.390625, 8985.390625 ]   (range: 17970.78125 )
%
% Effective contrast ratio (max range / min range): 615.052080197
%
% Known global minimum (WORK-space):
%   x* = [149289.0625000002;648.7890625000009;399039.0625;898.5390625]
%   f* = 0
%
% USAGE:
%   f = wood_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = wood_4D(n)   % Get bounds for dimension n (must be 4)
%   info = wood_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [10;10;10;10];
lb_work = [-1492890.625;-6487.890625;-3990390.625;-8985.390625];
ub_work = [1492890.625;6487.890625;3990390.625;8985.390625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625];
contrast_ratio = 615.0520801974833;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'wood';
    info.source_p = 63;
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
    info.x_global_min_orig = [1;1;1;1];
    info.x_global_min_work = [149289.0625000002;648.7890625000009;399039.0625;898.5390625];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Wood (4D): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = wood_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = wood_orig(x);
%
% Purpose:
%
%    Function wood is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [1 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = 100*(x(1,1)^2 - x(2,1))^2 + (1-x(1,1))^2 + 90*(x(3,1)^2-x(4,1))^2 +...
   (1-x(3,1))^2 + 10.1*((1-x(2,1))^2 + (1-x(4,1))^2) + 19.8*(1-x(2,1))*(1-x(4,1));
%
% End of wood.

