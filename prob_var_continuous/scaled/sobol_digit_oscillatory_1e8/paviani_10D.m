function varargout = paviani_10D(varargin)
%PAVIANI_10D  paviani 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [83353102.3595, 299908711.173]   (range: 216555608.814)
%   x2   ∈ [104310988.46, 375316253.682]   (range: 271005265.222)
%   x3   ∈ [10293346.0495, 37035983.7844]   (range: 26742637.7349)
%   x4   ∈ [17727.2226613, 63783.4507718]   (range: 46056.2281105)
%   x5   ∈ [58518125.8731, 210551199.824]   (range: 152033073.951)
%   x6   ∈ [12914.6346589, 46467.5138198]   (range: 33552.8791609)
%   x7   ∈ [3549.04805669, 12769.655819]   (range: 9220.60776228)
%   x8   ∈ [15207.5473692, 54717.5306317]   (range: 39509.9832625)
%   x9   ∈ [9660.29943085, 34758.2497813]   (range: 25097.9503504)
%   x10  ∈ [116071023.707, 417629459.963]   (range: 301558436.255)
%
% Effective contrast ratio (max range / min range): 32704.8329166
%
% USAGE:
%   f = paviani_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = paviani_10D(n)   % Get bounds for dimension n (must be 10)
%   info = paviani_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997];
lb_work = [83353102.35953546;104310988.4603721;10293346.04948725;17727.22266127964;58518125.87312391;12914.63465885557;3549.048056692566;15207.54736919898;9660.299430852967;116071023.7072489];
ub_work = [299908711.1733871;375316253.6822297;37035983.78435452;63783.45077182162;210551199.824453;46467.51381977133;12769.65581897525;54717.53063169512;34758.24978126543;417629459.9625587];
scale_factors = [41655723.31810868;52129429.51542838;5144100.974256499;8859.181739769938;29244440.71620386;6454.090284285641;1773.637209741413;7599.973697750616;4827.735847502734;58006508.59932479];
contrast_ratio = 32704.83291663792;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'paviani';
    info.source_p = 36;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
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
varargout{1} = paviani_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = paviani_orig(x);
%
% Purpose:
%
%    Function paviani is the function described in
%    Montaz et al. (2005).
%
%    dim = 10
%    Minimum global value = -45.778
%    Global minima = 9.351*ones(10,1)
%    Local minima
%    Search domain: 2.001 <= x <= 9.999
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
f = sum((log(x-2)).^2+(log(10-x)).^2,1)-prod(x,1)^0.2;
%
% End of paviani.

