function varargout = rosenbrock_2D(varargin)
%ROSENBROCK_2D  rosenbrock 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [-5.12       , 5.12        ]   (range: 10.24       )
%   x2   ∈ [-1536       , 1536        ]   (range: 3072        )
%
% Effective contrast ratio (max range / min range): 300
%
% Known global minimum (WORK-space):
%   x* = [1;300]
%   f* = 0
%
% USAGE:
%   f = rosenbrock_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = rosenbrock_2D(n)   % Get bounds for dimension n (must be 2)
%   info = rosenbrock_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-5.12;-5.12];
ub_orig = [5.12;5.12];
lb_work = [-5.12;-1536];
ub_work = [5.12;1536];
scale_factors = [1;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'rosenbrock';
    info.source_p = 41;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;1];
    info.x_global_min_work = [1;300];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Rosenbrock (n=2): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = rosenbrock_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = rosenbrock_orig(x);
%
% Purpose:
%
%    Function rosenbrock is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (equal to global minima)
%    Search domain: -1000 <= x <= 1000 (Brachetti et al. (1997))
%                   -2.048 <= x <= 2.048 (Storn and Price (1997))
%                   -5.12 <= x <= 5.12 (Huyer and Neumaier (2008))
%    Cases considered: n = 2 Huyer and Neumaier (2008)
%                            Brachetti et al. (1997)
%                      n = 10 Garcia-Palomares  et al. (2006)
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
xaux = x(2:n);
f    = sum((x(1:n-1) - 1).^2 + 100.*(x(1:n-1).^2 - xaux).^2,1);
%
% End of rosenbrock.

