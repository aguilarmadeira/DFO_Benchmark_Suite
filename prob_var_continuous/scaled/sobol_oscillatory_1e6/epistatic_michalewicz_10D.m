function varargout = epistatic_michalewicz_10D(varargin)
%EPISTATIC_MICHALEWICZ_10D  epistatic_michalewicz 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 304853.524307]   (range: 304853.524307)
%   x2   ∈ [0           , 1324.85011911]   (range: 1324.85011911)
%   x3   ∈ [0           , 814851.821709]   (range: 814851.821709)
%   x4   ∈ [0           , 1834.84841651]   (range: 1834.84841651)
%   x5   ∈ [0           , 177353.949957]   (range: 177353.949957)
%   x6   ∈ [0           , 1197.35054476]   (range: 1197.35054476)
%   x7   ∈ [0           , 687352.247359]   (range: 687352.247359)
%   x8   ∈ [0           , 1707.34884216]   (range: 1707.34884216)
%   x9   ∈ [0           , 432353.098658]   (range: 432353.098658)
%   x10  ∈ [0           , 1452.34969346]   (range: 1452.34969346)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   see info.x_global_min_work (not stored as a representative here)
%   f* = -9.66
%
% USAGE:
%   f = epistatic_michalewicz_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = epistatic_michalewicz_10D(n)   % Get bounds for dimension n (must be 10)
%   info = epistatic_michalewicz_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366;2.042035224833366];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [304853.5243073499;1324.850119111616;814851.8217094829;1834.848416513749;177353.9499568166;1197.350544761083;687352.2473589496;1707.348842163216;432353.0986578831;1452.349693462149];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'epistatic_michalewicz';
    info.source_p = 12;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -9.66;
    info.x_global_min_orig = [];
    info.x_global_min_work = [];
    info.global_min_note = 'Global minimizer is known but infoG.xstar_orig is missing/empty.';
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

