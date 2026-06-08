function varargout = wood_4D(varargin)
%WOOD_4D  wood 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-62381.5860456, 18714.4758137]   (range: 81096.0618592)
%   x2   ∈ [-11225.7773449, 3367.73320348]   (range: 14593.5105484)
%   x3   ∈ [-88755.732784, 26626.7198352]   (range: 115382.452619)
%   x4   ∈ [-37577.65541, 11273.296623]   (range: 48850.952033)
%
% Effective contrast ratio (max range / min range): 7.90642198369
%
% Known global minimum (WORK-space):
%   x* = [6238.158604555167;1122.577734492912;8875.573278401091;3757.765541003631]
%   f* = 0
%
% USAGE:
%   f = wood_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = wood_4D(n)   % Get bounds for dimension n (must be 4)
%   info = wood_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [3;3;3;3];
lb_work = [-62381.5860455517;-11225.77734492911;-88755.73278401091;-37577.65541003633];
ub_work = [18714.4758136655;3367.733203478734;26626.71983520327;11273.2966230109];
scale_factors = [6238.158604555169;1122.577734492911;8875.573278401091;3757.765541003632];
contrast_ratio = 7.90642198369483;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'wood';
    info.source_p = 63;
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
    info.x_global_min_orig = [1;1;1;1];
    info.x_global_min_work = [6238.158604555167;1122.577734492912;8875.573278401091;3757.765541003631];
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

