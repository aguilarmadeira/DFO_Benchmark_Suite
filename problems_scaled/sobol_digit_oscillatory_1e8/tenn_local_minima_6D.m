function varargout = tenn_local_minima_6D(varargin)
%TENN_LOCAL_MINIMA_6D  tenn_local_minima 6D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-749666069.107, 749666069.107]   (range: 1499332138.21)
%   x2   ∈ [-128216649.963, 128216649.963]   (range: 256433299.927)
%   x3   ∈ [-941322123.373, 941322123.373]   (range: 1882644246.75)
%   x4   ∈ [-437059675.582, 437059675.582]   (range: 874119351.164)
%   x5   ∈ [-584772542.52, 584772542.52]   (range: 1169545085.04)
%   x6   ∈ [-4143.95938717, 4143.95938717]   (range: 8287.91877433)
%
% Effective contrast ratio (max range / min range): 227155.248261
%
% Known global minimum (WORK-space):
%   x* = [74966606.91070354;12821664.99633081;94132212.33733416;43705967.55818456;58477254.25202155;414.3959387165733]
%   f* = 0
%
% USAGE:
%   f = tenn_local_minima_6D(x)          % Evaluate function at point x (6D vector)
%   [lb, ub] = tenn_local_minima_6D(n)   % Get bounds for dimension n (must be 6)
%   info = tenn_local_minima_6D()        % Get complete problem information

nloc = 6;
lb_orig = [-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10];
lb_work = [-749666069.1070352;-128216649.963308;-941322123.3733406;-437059675.5818454;-584772542.5202154;-4143.959387165729];
ub_work = [749666069.1070352;128216649.963308;941322123.3733406;437059675.5818454;584772542.5202154;4143.959387165729];
scale_factors = [74966606.91070351;12821664.9963308;94132212.33733407;43705967.55818454;58477254.25202154;414.3959387165729];
contrast_ratio = 227155.2482605676;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'tenn_local_minima';
    info.source_p = 59;
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
    info.x_global_min_orig = [1;1;1;1;1;1];
    info.x_global_min_work = [74966606.91070354;12821664.99633081;94132212.33733416;43705967.55818456;58477254.25202155;414.3959387165733];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Ten Local Minima (n=6): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 6D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 6 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = tenn_local_minima_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = tenn_local_minima_orig(x);
%
% Purpose:
%
%    Function tenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (10^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
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
f   = pi/dim*sum(((x(1:dim-1)-1).^2).*(1+10*sin(pi.*x(2:dim)).^2),1);
f   = f + pi/dim*(10*sin(pi*x(1))^2+(x(dim)-1)^2);
%
% End of tenn_local_minima.

