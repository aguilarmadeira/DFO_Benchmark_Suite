function varargout = becker_lago_2D(varargin)
%BECKER_LAGO_2D  becker_lago 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-5005000    , 5005000     ]   (range: 10010000    )
%   x2   ∈ [-2507.5     , 2507.5      ]   (range: 5015        )
%
% Effective contrast ratio (max range / min range): 1996.01196411
%
% Known global minimum (WORK-space):
%   multiple minimizers; see info.x_global_min_work
%   f* = 0
%
% USAGE:
%   f = becker_lago_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = becker_lago_2D(n)   % Get bounds for dimension n (must be 2)
%   info = becker_lago_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-10;-10];
ub_orig = [10;10];
lb_work = [-5005000;-2507.5];
ub_work = [5005000;2507.5];
scale_factors = [500500;250.75];
contrast_ratio = 1996.011964107677;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'becker_lago';
    info.source_p = 3;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [-5 -5 5 5;-5 5 -5 5];
    info.x_global_min_work = [-2502500 -2502500 2502500 2502500;-1253.75 1253.75 -1253.75 1253.75];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Becker-Lago (2D): 4 global minima at (+-5,+-5), f*=0. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 2D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 2 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = becker_lago_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = becker_lago_orig(x);
%
% Purpose:
%
%    Function becker_lago is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [-5 -5]', [-5 5]', [5 -5]', [5 5]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
f = (abs(x(1))-5)^2+(abs(x(2))-5)^2;
%
% End of becker_lago.

