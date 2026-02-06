function varargout = exponencial_2D(varargin)
%EXPONENCIAL_2D  exponencial 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-500500     , 500500      ]   (range: 1001000     )
%   x2   ∈ [-250.75     , 250.75      ]   (range: 501.5       )
%
% Effective contrast ratio (max range / min range): 1996.01196411
%
% Known global minimum (WORK-space):
%   x* = [0;0]
%   f* = -1
%
% USAGE:
%   f = exponencial_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = exponencial_2D(n)   % Get bounds for dimension n (must be 2)
%   info = exponencial_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-1;-1];
ub_orig = [1;1];
lb_work = [-500500;-250.75];
ub_work = [500500;250.75];
scale_factors = [500500;250.75];
contrast_ratio = 1996.011964107677;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'exponencial';
    info.source_p = 13;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 1;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1;
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Exponencial (n=2): x*=0, f*=-1. Ref: Brachetti et al. (1997).';
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
varargout{1} = exponencial_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = exponencial_orig(x);
%
% Purpose:
%
%    Function exponencial is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = -1
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -1 <= x <= 1 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
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
f = -exp(-0.5*sum(x.^2,1)); 
%
% End of exponencial.

