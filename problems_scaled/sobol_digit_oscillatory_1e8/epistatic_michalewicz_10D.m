function varargout = epistatic_michalewicz_10D(varargin)
%EPISTATIC_MICHALEWICZ_10D  epistatic_michalewicz 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 31303.5312806]   (range: 31303.5312806)
%   x2   ∈ [0           , 6477.1733931]   (range: 6477.1733931)
%   x3   ∈ [0           , 19590.920395]   (range: 19590.920395)
%   x4   ∈ [0           , 103672718.424]   (range: 103672718.424)
%   x5   ∈ [0           , 267471996.721]   (range: 267471996.721)
%   x6   ∈ [0           , 22294008.527]   (range: 22294008.527)
%   x7   ∈ [0           , 227381095.284]   (range: 227381095.284)
%   x8   ∈ [0           , 14061.00295 ]   (range: 14061.00295 )
%   x9   ∈ [0           , 28723.6761337]   (range: 28723.6761337)
%   x10  ∈ [0           , 4156.02804171]   (range: 4156.02804171)
%
% Effective contrast ratio (max range / min range): 64357.6015456
%
% Known global minimum (WORK-space):
%   see info.x_global_min_work (not stored as a representative here)
%   f* = -9.66
%
% USAGE:
%   f = epistatic_michalewicz_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = epistatic_michalewicz_10D(n)   % Get bounds for dimension n (must be 10)
%   info = epistatic_michalewicz_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793;3.141592653589793];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [31303.531280606;6477.173393097561;19590.92039496387;103672718.4242032;267471996.7205379;22294008.5270444;227381095.2840327;14061.00294996852;28723.67613373026;4156.028041705326];
scale_factors = [9964.223479080429;2061.748325549562;6235.983641156654;33000051.20197229;85138980.83346567;7096403.316823962;72377650.56020612;4475.756248634426;9143.030080907743;1322.904812931865];
contrast_ratio = 64357.60154563039;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'epistatic_michalewicz';
    info.source_p = 12;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -9.66;
    info.x_global_min_orig = [];
    info.x_global_min_work = [];
    info.global_min_note = 'Global minimizer is known but infoG.xstar_orig is missing/empty.';
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
varargout{1} = epistatic_michalewicz_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = epistatic_michalewicz_orig(x);
%
% Purpose:
%
%    Function epistatic_michalewicz is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -4.687 if n = 5
%                          -9.66 if n = 10
%    Global minima
%    Local minima (several)
%    Search domain: 0 <= x <= pi (Huyer and Neumaier (1999))
%    Cases considered: n = 5 Huyer and Neumaier (1999), ICEO
%                      n = 10 Huyer and Neumaier (1999), ICEO
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
n    = size(x,1);
even = 0;
for i=1:n-1
    if even == 0
       y(i,1) = x(i)*cos(pi/6)-x(i+1)*sin(pi/6);
       even = 1;
    else
       y(i,1) = x(i)*sin(pi/6)+x(i+1)*cos(pi/6);
       even = 0;
    end
end
y = [y;x(n)];
f = -sum(sin(y).*(sin([1:n]'.*(y.^2)./pi)).^20,1);
%
% End of epistatic_michalewicz.

