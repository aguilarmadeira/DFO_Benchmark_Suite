function varargout = epistatic_michalewicz_5D(varargin)
%EPISTATIC_MICHALEWICZ_5D  epistatic_michalewicz 5D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 304853.524307]   (range: 304853.524307)
%   x2   ∈ [0           , 1324.85011911]   (range: 1324.85011911)
%   x3   ∈ [0           , 814851.821709]   (range: 814851.821709)
%   x4   ∈ [0           , 1834.84841651]   (range: 1834.84841651)
%   x5   ∈ [0           , 177353.949957]   (range: 177353.949957)
%
% Effective contrast ratio (max range / min range): 615.052080197
%
% Known global minimum (WORK-space):
%   see info.x_global_min_work (not stored as a representative here)
%   f* = -4.687
%
% USAGE:
%   f = epistatic_michalewicz_5D(x)          % Evaluate function at point x (5D vector)
%   [lb, ub] = epistatic_michalewicz_5D(n)   % Get bounds for dimension n (must be 5)
%   info = epistatic_michalewicz_5D()        % Get complete problem information

nloc = 5;
lb_orig = [0;0;0;0;0];
ub_orig = [2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366];
lb_work = [0;0;0;0;0];
ub_work = [304853.5243073499;1324.850119111616;814851.8217094829;1834.848416513749;177353.9499568166];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625];
contrast_ratio = 615.0520801974833;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'epistatic_michalewicz';
    info.source_p = 11;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -4.687;
    info.x_global_min_orig = [];
    info.x_global_min_work = [];
    info.global_min_note = 'Global minimizer is known but infoG.xstar_orig is missing/empty.';
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
varargout{1} = epistatic_michalewicz_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = epistatic_michalewicz_orig(x);
%
% Purpose:
%
%    Function epistatic_michalewicz is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -4.687 if n = 5
%                          -9.66 if n = 10
%    Global minima
%    Local minima (several)
%    Search domain: 0 <= x <= pi (Huyer and Neumaier (1999))
%    Cases considered: n = 5 Huyer and Neumaier (1999), ICEO
%                      n = 10 Huyer and Neumaier (1999), ICEO
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
even = 0;
for i=1:n-1
    if even == 0
       y(i,1) = x(i)*cos(pi/6)-x(i+1)*sin(pi/6);
       even = 1;
    else
       y(i,1) = x(i)*sin(pi/6)+x(i+1)*cos(pi/6);
       even = 0;
    end
end
y = [y;x(n)];
f = -sum(sin(y).*(sin([1:n]'.*(y.^2)./pi)).^20,1);
%
% End of epistatic_michalewicz.

