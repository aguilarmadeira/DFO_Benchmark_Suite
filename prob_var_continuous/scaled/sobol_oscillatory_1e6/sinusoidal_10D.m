function varargout = sinusoidal_10D(varargin)
%SINUSOIDAL_10D  sinusoidal 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 17466820.3125]   (range: 17466820.3125)
%   x2   ∈ [0           , 75908.3203125]   (range: 75908.3203125)
%   x3   ∈ [0           , 46687570.3125]   (range: 46687570.3125)
%   x4   ∈ [0           , 105129.070312]   (range: 105129.070312)
%   x5   ∈ [0           , 10161632.8125]   (range: 10161632.8125)
%   x6   ∈ [0           , 68603.1328125]   (range: 68603.1328125)
%   x7   ∈ [0           , 39382382.8125]   (range: 39382382.8125)
%   x8   ∈ [0           , 97823.8828125]   (range: 97823.8828125)
%   x9   ∈ [0           , 24772007.8125]   (range: 24772007.8125)
%   x10  ∈ [0           , 83213.5078125]   (range: 83213.5078125)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% USAGE:
%   f = sinusoidal_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = sinusoidal_10D(n)   % Get bounds for dimension n (must be 10)
%   info = sinusoidal_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [117;117;117;117;117;117;117;117;117;117];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [17466820.3125;75908.3203125;46687570.3125;105129.0703125;10161632.8125;68603.1328125;39382382.8125;97823.8828125;24772007.8125;83213.5078125];
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
    info.problem = 'sinusoidal';
    info.source_p = 53;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
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

