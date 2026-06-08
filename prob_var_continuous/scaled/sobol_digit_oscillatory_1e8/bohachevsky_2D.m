function varargout = bohachevsky_2D(varargin)
%BOHACHEVSKY_2D  bohachevsky 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-269424.927756, 80827.4783269]   (range: 350252.406083)
%   x2   ∈ [-740302954.633, 222090886.39]   (range: 962393841.024)
%
% Effective contrast ratio (max range / min range): 2747.71514573
%
% Known global minimum (WORK-space):
%   x* = [0;1.192092895507812e-07]
%   f* = 0
%
% USAGE:
%   f = bohachevsky_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = bohachevsky_2D(n)   % Get bounds for dimension n (must be 2)
%   info = bohachevsky_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-50;-50];
ub_orig = [15;15];
lb_work = [-269424.9277564717;-740302954.6334758];
ub_work = [80827.47832694151;222090886.3900427];
scale_factors = [5388.498555129433;14806059.09266951];
contrast_ratio = 2747.715145729286;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'bohachevsky';
    info.source_p = 4;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 50;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;1.192092895507812e-07];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Bohachevsky (2D): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = bohachevsky_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = bohachevsky_orig(x);
%
% Purpose:
%
%    Function bohachevsky is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -50 <= x <= 50 (Montaz et al. (2005))
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
f = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2))+0.3;
%
% End of bohachevsky.

