function varargout = cauchy_10D(varargin)
%CAUCHY_10D  cauchy 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [905.310927512, 7966.73616211]   (range: 7061.4252346)
%   x2   ∈ [13620.1910159, 119857.68094]   (range: 106237.489924)
%   x3   ∈ [9691.07295222, 85281.4419795]   (range: 75590.3690273)
%   x4   ∈ [170158636.762, 1497396003.51]   (range: 1327237366.75)
%   x5   ∈ [4932.17633802, 43403.1517746]   (range: 38470.9754366)
%   x6   ∈ [12097.9960648, 106462.36537]   (range: 94364.3693055)
%   x7   ∈ [5666.95384072, 49869.1937983]   (range: 44202.2399576)
%   x8   ∈ [185353123.124, 1631107483.49]   (range: 1445754360.36)
%   x9   ∈ [24235847.0062, 213275453.654]   (range: 189039606.648)
%   x10  ∈ [14627.9478148, 128725.94077]   (range: 114097.992955)
%
% Effective contrast ratio (max range / min range): 204739.739122
%
% USAGE:
%   f = cauchy_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = cauchy_10D(n)   % Get bounds for dimension n (must be 10)
%   info = cauchy_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2;2;2;2;2;2;2;2;2;2];
ub_orig = [17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6;17.6];
lb_work = [905.3109275122283;13620.19101594932;9691.072952217859;170158636.7623068;4932.176338020254;12097.99606480983;5666.953840715239;185353123.1235698;24235847.00617321;14627.94781476223];
ub_work = [7966.736162107612;119857.680940354;85281.44197951719;1497396003.5083;43403.15177457823;106462.3653703266;49869.19379829412;1631107483.487414;213275453.6543242;128725.9407699077];
scale_factors = [452.6554637561143;6810.095507974661;4845.53647610893;85079318.38115339;2466.088169010127;6048.998032404917;2833.47692035762;92676561.56178491;12117923.5030866;7313.973907381116];
contrast_ratio = 204739.7391224643;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'cauchy';
    info.source_p = 7;
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
varargout{1} = cauchy_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = cauchy_orig(x);
%
% Purpose:
%
%    Function cauchy is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value
%    Global minima
%    Local minima
%    Search domain: y(1) <= x <= y(dim) (Brachetti et al. (1997))
%    Cases considered: n = 4 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
%                      n = 25 Brachetti et al. (1997)
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
n = size(x,1);
if n == 4
    y = [3 7 12 17]';
else
    if n == 10
        y = [2 5 7 8 11 15 17 21 23 26]';
    else
        if n == 25
            y = [4.1 7.7 17.5 31.4 32.7 92.4 115.3 118.3 119 129.6 198.6];
            y = [y, 200.7 242.5 255 274.7 274.7 303.8 334.1 430 489.1 703.4 ];
            y = [y, 978 1656 1697.8 2745.6]';
        end
    end
end
f = -sum(log(pi)+log(1+(y-x).^2),1);
%
% End of cauchy.

