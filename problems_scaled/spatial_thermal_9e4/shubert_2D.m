function varargout = shubert_2D(varargin)
%SHUBERT_2D  shubert 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SPATIAL THERMAL HETEROGENEITY):
%
%   x1   ∈ [-10         , 10          ]   (range: 20          )
%   x2   ∈ [-3000       , 3000        ]   (range: 6000        )
%
% Effective contrast ratio (max range / min range): 300
%
% Known global minimum (WORK-space):
%   see info.x_global_min_work (not stored as a representative here)
%   f* = -186.7309
%
% USAGE:
%   f = shubert_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = shubert_2D(n)   % Get bounds for dimension n (must be 2)
%   info = shubert_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-10;-10];
ub_orig = [10;10];
lb_work = [-10;-3000];
ub_work = [10;3000];
scale_factors = [1;300];
contrast_ratio = 300;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'shubert';
    info.source_p = 52;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -186.7309;
    info.x_global_min_orig = [];
    info.x_global_min_work = [];
    info.global_min_note = 'Global minimizer is known but infoG.xstar_orig is missing/empty.';
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
varargout{1} = shubert_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = shubert_orig(x);
%
% Purpose:
%
%    Function shubert is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -186.7309
%    Global minima (18)
%    Local minima (760)
%    Search domain: -10 <= x <= 10 (Huyer and Neumaier (1999))
%                                  (Huyer and Neumaier (2008))
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
f = sum([1:5].*cos([2:6].*x(1)+[1:5]),2)*sum([1:5].*cos([2:6].*x(2)+[1:5]),2);
%
% End of shubert.

