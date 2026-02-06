function varargout = salomon_10D(varargin)
%SALOMON_10D  salomon 10D test problem (heterogeneous WORK-space wrapper).
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
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = salomon_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = salomon_10D(n)   % Get bounds for dimension n (must be 10)
%   info = salomon_10D()        % Get complete problem information

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
    info.problem = 'salomon';
    info.source_p = 43;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Salomon (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = salomon_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x);
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
%
% End of salomon.

