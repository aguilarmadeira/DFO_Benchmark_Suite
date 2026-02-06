function varargout = paviani_10D(varargin)
%PAVIANI_10D  paviani 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [2.001       , 9.999       ]   (range: 7.998       )
%   x2   ∈ [2.001       , 9.999       ]   (range: 7.998       )
%   x3   ∈ [2.001       , 9.999       ]   (range: 7.998       )
%   x4   ∈ [2.001       , 9.999       ]   (range: 7.998       )
%   x5   ∈ [2.001       , 9.999       ]   (range: 7.998       )
%   x6   ∈ [600.3       , 2999.7      ]   (range: 2399.4      )
%   x7   ∈ [600.3       , 2999.7      ]   (range: 2399.4      )
%   x8   ∈ [600.3       , 2999.7      ]   (range: 2399.4      )
%   x9   ∈ [600.3       , 2999.7      ]   (range: 2399.4      )
%   x10  ∈ [600.3       , 2999.7      ]   (range: 2399.4      )
%
% Effective contrast ratio (max range / min range): 300
%
% Known global minimum (WORK-space):
%   x* = [9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;2805.3;2805.3;2805.3;2805.3;2805.3]
%   f* = -45.778
%
% USAGE:
%   f = paviani_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = paviani_10D(n)   % Get bounds for dimension n (must be 10)
%   info = paviani_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001];
lb_work = [2.000999999999999;2.000999999999999;2.000999999999999;2.000999999999999;2.000999999999999;600.2999999999997;600.2999999999997;600.2999999999997;600.2999999999997;600.2999999999997];
ub_work = [9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;2999.7;2999.7;2999.7;2999.7;2999.7];
scale_factors = [1;1;1;1;1;300;300;300;300;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'paviani';
    info.source_p = 36;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -45.778;
    info.x_global_min_orig = [9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001];
    info.x_global_min_work = [9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;2805.3;2805.3;2805.3;2805.3;2805.3];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Paviani (10D): x*=9.351, f*=-45.778. Ref: Ali et al. (2005).';
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
varargout{1} = paviani_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = paviani_orig(x);
%
% Purpose:
%
%    Function paviani is the function described in
%    Montaz et al. (2005).
%
%    dim = 10
%    Minimum global value = -45.778
%    Global minima = 9.351*ones(10,1)
%    Local minima
%    Search domain: 2.001 <= x <= 9.999
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
f = sum((log(x-2)).^2+(log(10-x)).^2,1)-prod(x,1)^0.2;
%
% End of paviani.

