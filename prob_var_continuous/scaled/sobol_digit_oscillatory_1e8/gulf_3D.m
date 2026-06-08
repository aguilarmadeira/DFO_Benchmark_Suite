function varargout = gulf_3D(varargin)
%GULF_3D  gulf 3D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [775.092292628, 504081.272511]   (range: 503306.180218)
%   x2   ∈ [0           , 424407866.859]   (range: 424407866.859)
%   x3   ∈ [0           , 211345659.655]   (range: 211345659.655)
%
% Effective contrast ratio (max range / min range): 8389.89554645
%
% USAGE:
%   f = gulf_3D(x)          % Evaluate function at point x (3D vector)
%   [lb, ub] = gulf_3D(n)   % Get bounds for dimension n (must be 3)
%   info = gulf_3D()        % Get complete problem information

nloc = 3;
lb_orig = [0.1;0;0];
ub_orig = [65.035;16.64;3.25];
lb_work = [775.09229262796;0;0];
ub_work = [504081.2725106244;424407866.85862;211345659.6551653];
scale_factors = [7750.922926280072;25505280.4602536;65029433.74005087];
contrast_ratio = 8389.895546446967;

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
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
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

