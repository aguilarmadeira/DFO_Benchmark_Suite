function varargout = sinusoidal_10D(varargin)
%SINUSOIDAL_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   sinusoidal (source instance p=53)
% Dimension: n = 10
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 118572.6386671073
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = sinusoidal_orig(x_orig)

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [180;180;180;180;180;180;180;180;180;180];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [1519391.588670009;831536.7242391757;10086438034.66121;234637.067022919;1644180.348928585;506370.5113301522;13408405632.61487;116878.8536247006;13858634078.66735;7538120224.458486];
scale_factors = [8441.064381500051;4619.64846799542;56035766.85922894;1303.539261238439;9134.335271825474;2813.169507389734;74491142.40341596;649.32696458167;76992411.54815196;41878445.69143603];
contrast_ratio = 118572.6386671073;

if nargin == 0
    info.name = mfilename;
    info.problem = 'sinusoidal';
    info.source_p = 53;
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
varargout{1} = sinusoidal_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = sinusoidal_orig(x);
%
% Purpose:
%
%    Function sinusoidal is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -(A+1)
%    Global minima = (90+z)*ones(n,1)
%    Local minima (several)
%    Search domain: 0 <= x <= 180 (Montaz et al. (2005))
%    Cases considered: n = 10, 20
%                      A = 2.5;B = 5;z = 30;Montaz et al. (2005)
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
A = 2.5;
B = 5;
z = 30;
f = -(A*prod(sind(x-z),1)+prod(sind(B.*(x-z)),1));
%
% End of sinusoidal.

