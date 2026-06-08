function varargout = schaffer2_2D(varargin)
%SCHAFFER2_2D  schaffer2 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-576848.201071, 173054.460321]   (range: 749902.661393)
%   x2   ∈ [-2149564211.42, 644869263.427]   (range: 2794433474.85)
%
% Effective contrast ratio (max range / min range): 3726.39492926
%
% Known global minimum (WORK-space):
%   x* = [0;0]
%   f* = 0
%
% USAGE:
%   f = schaffer2_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = schaffer2_2D(n)   % Get bounds for dimension n (must be 2)
%   info = schaffer2_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-100;-100];
ub_orig = [30;30];
lb_work = [-576848.2010714116;-2149564211.424759];
ub_work = [173054.4603214235;644869263.4274276];
scale_factors = [5768.482010714116;21495642.11424759];
contrast_ratio = 3726.39492925913;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'schaffer2';
    info.source_p = 45;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Schaffer2 (2D): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = schaffer2_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schaffer2_orig(x);
%
% Purpose:
%
%    Function schaffer2 is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
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
f = (x(1)^2+x(2)^2)^0.25*(sin(50*(x(1)^2+x(2)^2)^0.1)^2+1);
%
% End of schaffer2.

