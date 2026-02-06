function varargout = neumaier3_10D(varargin)
%NEUMAIER3_10D  neumaier3 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (PROGRESSIVE HETEROGENEITY):
%
%   x1   ∈ [-100        , 100         ]   (range: 200         )
%   x2   ∈ [-1000       , 1000        ]   (range: 2000        )
%   x3   ∈ [-10000      , 10000       ]   (range: 20000       )
%   x4   ∈ [-100000     , 100000      ]   (range: 200000      )
%   x5   ∈ [-10         , 10          ]   (range: 20          )
%   x6   ∈ [-1          , 1           ]   (range: 2           )
%   x7   ∈ [-0.1        , 0.1         ]   (range: 0.2         )
%   x8   ∈ [-100        , 100         ]   (range: 200         )
%   x9   ∈ [-1000       , 1000        ]   (range: 2000        )
%   x10  ∈ [-10000      , 10000       ]   (range: 20000       )
%
% Effective contrast ratio (max range / min range): 1000000
%
% Known global minimum (WORK-space):
%   x* = [10.00000000000001;180;2400;28000;3;0.3;0.028;24;180;1000]
%   f* = -210
%
% USAGE:
%   f = neumaier3_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = neumaier3_10D(n)   % Get bounds for dimension n (must be 10)
%   info = neumaier3_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [100;100;100;100;100;100;100;100;100;100];
lb_work = [-100;-1000;-10000;-100000;-10;-1;-0.1;-100;-1000;-10000];
ub_work = [100;1000;10000;100000;10;1;0.1;100;1000;10000];
scale_factors = [1;10;100;1000;0.1;0.01;0.001;1;10;100];
contrast_ratio = 1000000;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'neumaier3';
    info.source_p = 35;
    info.dimension = nloc;
    info.strategy = 'progressive';
    info.kappa = 1000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -210;
    info.x_global_min_orig = [10;18;24;28;30;30;28;24;18;10];
    info.x_global_min_work = [10.00000000000001;180;2400;28000;3;0.3;0.028;24;180;1000];
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

