function varargout = salomon_5D(varargin)
%SALOMON_5D  salomon 5D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-1418622372.4, 425586711.72]   (range: 1844209084.12)
%   x2   ∈ [-795901.741245, 238770.522374]   (range: 1034672.26362)
%   x3   ∈ [-3246476264.36, 973942879.307]   (range: 4220419143.66)
%   x4   ∈ [-737853.242638, 221355.972791]   (range: 959209.21543)
%   x5   ∈ [-44151.5123827, 13245.4537148]   (range: 57396.9660975)
%
% Effective contrast ratio (max range / min range): 73530.3523969
%
% Known global minimum (WORK-space):
%   x* = [0;0;4.76837158203125e-07;1.164153218269348e-10;0]
%   f* = 0
%
% USAGE:
%   f = salomon_5D(x)          % Evaluate function at point x (5D vector)
%   [lb, ub] = salomon_5D(n)   % Get bounds for dimension n (must be 5)
%   info = salomon_5D()        % Get complete problem information

nloc = 5;
lb_orig = [-100;-100;-100;-100;-100];
ub_orig = [30;30;30;30;30];
lb_work = [-1418622372.398416;-795901.7412454954;-3246476264.355166;-737853.2426383159;-44151.51238271069];
ub_work = [425586711.7195248;238770.5223736486;973942879.3065498;221355.9727914948;13245.45371481321];
scale_factors = [14186223.72398416;7959.017412454954;32464762.64355166;7378.532426383159;441.5151238271069];
contrast_ratio = 73530.35239686276;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'salomon';
    info.source_p = 42;
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
    info.x_global_min_orig = [0;0;0;0;0];
    info.x_global_min_work = [0;0;4.76837158203125e-07;1.164153218269348e-10;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Salomon (n=5): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = salomon_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x);
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
%
% End of salomon.

