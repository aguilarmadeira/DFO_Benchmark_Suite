function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  schwefel 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-74644531.25, 22393359.375]   (range: 97037890.625)
%   x2   ∈ [-324394.53125, 97318.359375]   (range: 421712.890625)
%   x3   ∈ [-199519531.25, 59855859.375]   (range: 259375390.625)
%   x4   ∈ [-449269.53125, 134780.859375]   (range: 584050.390625)
%   x5   ∈ [-43425781.25, 13027734.375]   (range: 56453515.625)
%   x6   ∈ [-293175.78125, 87952.734375]   (range: 381128.515625)
%   x7   ∈ [-168300781.25, 50490234.375]   (range: 218791015.625)
%   x8   ∈ [-418050.78125, 125415.234375]   (range: 543466.015625)
%   x9   ∈ [-105863281.25, 31758984.375]   (range: 137622265.625)
%   x10  ∈ [-355613.28125, 106683.984375]   (range: 462297.265625)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% USAGE:
%   f = schwefel_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = schwefel_10D(n)   % Get bounds for dimension n (must be 10)
%   info = schwefel_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [150;150;150;150;150;150;150;150;150;150];
lb_work = [-74644531.25;-324394.53125;-199519531.25;-449269.53125;-43425781.25;-293175.78125;-168300781.25;-418050.78125;-105863281.25;-355613.28125];
ub_work = [22393359.375;97318.359375;59855859.375;134780.859375;13027734.375;87952.734375;50490234.375;125415.234375;31758984.375;106683.984375];
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
    info.problem = 'schwefel';
    info.source_p = 46;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 500;  % Original bound(p) from test setup
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
varargout{1} = schwefel_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schwefel_orig(x);
%
% Purpose:
%
%    Function schwefel is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -418.9829*n
%    Global minima = 420.97*ones(n,1)
%    Local minima 
%    Search domain: -500 <= x <= 500 (Montaz et al. (2005))
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
f = -sum(x.*sin(sqrt(abs(x))),1);
%
% End of schwefel.

