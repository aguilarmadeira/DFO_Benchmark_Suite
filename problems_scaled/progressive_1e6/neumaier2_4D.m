function varargout = neumaier2_4D(varargin)
%NEUMAIER2_4D  neumaier2 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (PROGRESSIVE HETEROGENEITY):
%
%   x1   ∈ [0           , 4           ]   (range: 4           )
%   x2   ∈ [0           , 40          ]   (range: 40          )
%   x3   ∈ [0           , 400         ]   (range: 400         )
%   x4   ∈ [0           , 4000        ]   (range: 4000        )
%
% Effective contrast ratio (max range / min range): 1000
%
% Known global minimum (WORK-space):
%   x* = [1;20;200;3000]
%   f* = 0
%
% USAGE:
%   f = neumaier2_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = neumaier2_4D(n)   % Get bounds for dimension n (must be 4)
%   info = neumaier2_4D()        % Get complete problem information

nloc = 4;
lb_orig = [0;0;0;0];
ub_orig = [4;4;4;4];
lb_work = [0;0;0;0];
ub_work = [4;40;400;4000];
scale_factors = [1;10;100;1000];
contrast_ratio = 1000;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'neumaier2';
    info.source_p = 34;
    info.dimension = nloc;
    info.strategy = 'progressive';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;2;2;3];
    info.x_global_min_work = [1;20;200;3000];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Neumaier2 (4D): x*=(1,2,2,3), f*=0. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 4D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 4 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = neumaier2_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = neumaier2_orig(x);
%
% Purpose:
%
%    Function neumaier2 is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = [1 2 2 3]'
%    Local minima
%    Search domain: 0 <= x <= n   (Montaz et al. (2005))
%    Cases considered: n = 4, b = [8 18 44 114]' Montaz et al. (2005)
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
b   = [8 18 44 114]';
aux = repmat([1:4]',1,4);
f   = sum((b - sum(repmat(x',4,1).^aux,2)).^2,1);
%
% End of neumaier2.

