function varargout = sinusoidal_10D(varargin)
%SINUSOIDAL_10D  sinusoidal 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 17347043853.8]   (range: 17347043853.8)
%   x2   ∈ [0           , 605252.546027]   (range: 605252.546027)
%   x3   ∈ [0           , 9009853761.66]   (range: 9009853761.66)
%   x4   ∈ [0           , 2288591501.08]   (range: 2288591501.08)
%   x5   ∈ [0           , 15430712295.3]   (range: 15430712295.3)
%   x6   ∈ [0           , 870725.298892]   (range: 870725.298892)
%   x7   ∈ [0           , 12226819344.4]   (range: 12226819344.4)
%   x8   ∈ [0           , 934497996.665]   (range: 934497996.665)
%   x9   ∈ [0           , 1597037.88435]   (range: 1597037.88435)
%   x10  ∈ [0           , 474960.094221]   (range: 474960.094221)
%
% Effective contrast ratio (max range / min range): 36523.1607136
%
% USAGE:
%   f = sinusoidal_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = sinusoidal_10D(n)   % Get bounds for dimension n (must be 10)
%   info = sinusoidal_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [180;180;180;180;180;180;180;180;180;180];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [17347043853.79442;605252.5460272937;9009853761.660627;2288591501.079745;15430712295.33369;870725.298891934;12226819344.35275;934497996.6647424;1597037.884350121;474960.0942212867];
scale_factors = [96372465.85441345;3362.514144596076;50054743.12033682;12714397.22822081;85726179.41852051;4837.362771621855;67926774.13529308;5191655.537026347;8872.432690834003;2638.667190118259];
contrast_ratio = 36523.16071360793;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'sinusoidal';
    info.source_p = 53;
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
varargout{1} = sinusoidal_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = sinusoidal_orig(x);
%
% Purpose:
%
%    Function sinusoidal is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -(A+1)
%    Global minima = (90+z)*ones(n,1)
%    Local minima (several)
%    Search domain: 0 <= x <= 180 (Montaz et al. (2005))
%    Cases considered: n = 10, 20
%                      A = 2.5;B = 5;z = 30;Montaz et al. (2005)
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
A = 2.5;
B = 5;
z = 30;
f = -(A*prod(sind(x-z),1)+prod(sind(B.*(x-z)),1));
%
% End of sinusoidal.

