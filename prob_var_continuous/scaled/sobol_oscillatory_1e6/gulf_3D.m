function varargout = gulf_3D(varargin)
%GULF_3D  gulf 3D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [14928.90625 , 9709014.17969]   (range: 9694085.27344)
%   x2   ∈ [0           , 10795.85    ]   (range: 10795.85    )
%   x3   ∈ [0           , 1296876.95312]   (range: 1296876.95312)
%
% Effective contrast ratio (max range / min range): 615.052080197
%
% USAGE:
%   f = gulf_3D(x)          % Evaluate function at point x (3D vector)
%   [lb, ub] = gulf_3D(n)   % Get bounds for dimension n (must be 3)
%   info = gulf_3D()        % Get complete problem information

nloc = 3;
lb_orig = [0.1;0;0];
ub_orig = [65.035;16.64;3.25];
lb_work = [14928.90624999907;0;0];
ub_work = [9709014.1796875;10795.85;1296876.953125];
scale_factors = [149289.0625;648.7890625;399039.0625];
contrast_ratio = 615.0520801974833;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'gulf';
    info.source_p = 24;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = false;
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 3D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 3 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = gulf_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = gulf_orig(x);
%
% Purpose:
%
%    Function gulf is the function described in
%    Montaz et al. (2005).
%
%    dim = 3
%    Minimum global value = 0
%    Global minima = [50 25 1.5]'
%    Local minima 
%    Search domain: 0.1 <= x(1) <= 100 
%                   0 <= x(2) <= 25.6
%                   0 <= x(3) <= 5    (Montaz et al. (2005))
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
u = 25+(-50*log(0.01.*[1:99])).^(2/3);
f = sum((exp(-(u-x(2)).^x(3)./x(1))-0.01.*[1:99]).^2,2);
%
% End of gulf.

