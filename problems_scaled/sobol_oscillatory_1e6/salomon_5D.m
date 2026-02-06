function varargout = salomon_5D(varargin)
%SALOMON_5D  salomon 5D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-14928906.25, 14928906.25 ]   (range: 29857812.5  )
%   x2   ∈ [-64878.90625, 64878.90625 ]   (range: 129757.8125 )
%   x3   ∈ [-39903906.25, 39903906.25 ]   (range: 79807812.5  )
%   x4   ∈ [-89853.90625, 89853.90625 ]   (range: 179707.8125 )
%   x5   ∈ [-8685156.25 , 8685156.25  ]   (range: 17370312.5  )
%
% Effective contrast ratio (max range / min range): 615.052080197
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = salomon_5D(x)          % Evaluate function at point x (5D vector)
%   [lb, ub] = salomon_5D(n)   % Get bounds for dimension n (must be 5)
%   info = salomon_5D()        % Get complete problem information

nloc = 5;
lb_orig = [-100;-100;-100;-100;-100];
ub_orig = [100;100;100;100;100];
lb_work = [-14928906.25;-64878.90625;-39903906.25;-89853.90625;-8685156.25];
ub_work = [14928906.25;64878.90625;39903906.25;89853.90625;8685156.25];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625];
contrast_ratio = 615.0520801974833;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'salomon';
    info.source_p = 42;
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
    info.x_global_min_orig = [0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Salomon (n=5): x*=0, f*=0. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 5D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 5 components.');
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

