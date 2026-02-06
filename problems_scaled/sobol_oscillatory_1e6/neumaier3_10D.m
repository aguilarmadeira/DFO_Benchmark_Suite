function varargout = neumaier3_10D(varargin)
%NEUMAIER3_10D  neumaier3 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-14928906.25, 14928906.25 ]   (range: 29857812.5  )
%   x2   ∈ [-64878.90625, 64878.90625 ]   (range: 129757.8125 )
%   x3   ∈ [-39903906.25, 39903906.25 ]   (range: 79807812.5  )
%   x4   ∈ [-89853.90625, 89853.90625 ]   (range: 179707.8125 )
%   x5   ∈ [-8685156.25 , 8685156.25  ]   (range: 17370312.5  )
%   x6   ∈ [-58635.15625, 58635.15625 ]   (range: 117270.3125 )
%   x7   ∈ [-33660156.25, 33660156.25 ]   (range: 67320312.5  )
%   x8   ∈ [-83610.15625, 83610.15625 ]   (range: 167220.3125 )
%   x9   ∈ [-21172656.25, 21172656.25 ]   (range: 42345312.5  )
%   x10  ∈ [-71122.65625, 71122.65625 ]   (range: 142245.3125 )
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   x* = [1492890.625000002;11678.203125;9576937.5;25159.09375;2605546.875;17590.546875;9424843.75;20066.4375;3811078.125;7112.265625]
%   f* = -210
%
% USAGE:
%   f = neumaier3_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = neumaier3_10D(n)   % Get bounds for dimension n (must be 10)
%   info = neumaier3_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [100;100;100;100;100;100;100;100;100;100];
lb_work = [-14928906.25;-64878.90625;-39903906.25;-89853.90625;-8685156.25;-58635.15625;-33660156.25;-83610.15625;-21172656.25;-71122.65625];
ub_work = [14928906.25;64878.90625;39903906.25;89853.90625;8685156.25;58635.15625;33660156.25;83610.15625;21172656.25;71122.65625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'neumaier3';
    info.source_p = 35;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -210;
    info.x_global_min_orig = [10;18;24;28;30;30;28;24;18;10];
    info.x_global_min_work = [1492890.625000002;11678.203125;9576937.5;25159.09375;2605546.875;17590.546875;9424843.75;20066.4375;3811078.125;7112.265625];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Neumaier3 (10D): x*_i=i*(n+1-i), f*=-210. Ref: Ali et al. (2005).';
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
varargout{1} = neumaier3_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = neumaier3_orig(x);
%
% Purpose:
%
%    Function neumaier3 is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -n*(n+4)*(n-1)/6
%    Global minima = [1:n]'.*(n+1-[1:n]')
%    Local minima (several)
%    Search domain: -n^2 <= x <= n^2   (Montaz et al. (2005))
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
dim  = size(x,1);
xaux = x(1:dim-1,1);
f    = sum((x-1).^2,1)-sum(x(2:dim,1).*xaux,1);
%
% End of neumaier3.

