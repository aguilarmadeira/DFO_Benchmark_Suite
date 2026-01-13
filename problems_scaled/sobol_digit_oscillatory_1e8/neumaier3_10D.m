function varargout = neumaier3_10D(varargin)
%NEUMAIER3_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   neumaier3 (source instance p=35)
% Dimension: n = 10
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 100
% Effective contrast: 40482.73457099601
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = neumaier3_orig(x_orig)

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [100;100;100;100;100;100;100;100;100;100];
lb_work = [-824545.8138585094;-2514562799.046238;-6167105700.75314;-569940882.9282978;-8948659829.358858;-462395.0722741366;-6714089110.416155;-221048.7983133964;-7801846072.519416;-331952.8842256456];
ub_work = [824545.8138585094;2514562799.046238;6167105700.75314;569940882.9282978;8948659829.358858;462395.0722741366;6714089110.416155;221048.7983133964;7801846072.519416;331952.8842256456];
scale_factors = [8245.458138585094;25145627.99046238;61671057.0075314;5699408.829282978;89486598.29358858;4623.950722741366;67140891.10416155;2210.487983133964;78018460.72519416;3319.528842256456];
contrast_ratio = 40482.73457099601;

if nargin == 0
    info.name = mfilename;
    info.problem = 'neumaier3';
    info.source_p = 35;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -210;
    info.x_global_min_orig = [10;18;24;28;30;30;28;24;18;10];
    info.x_global_min_work = [82454.58138585102;452621303.8283229;1480105368.180754;159583447.2199234;2684597948.807657;138718.521682241;1879944950.916523;53051.71159521514;1404332293.053494;33195.28842256457];
    info.global_min_note = 'Neumaier3 (10D): x*_i=i*(n+1-i), f*=-210. Ref: Ali et al. (2005).';
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
varargout{1} = neumaier3_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = neumaier3_orig(x);
%
% Purpose:
%
%    Function neumaier3 is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -n*(n+4)*(n-1)/6
%    Global minima = [1:n]'.*(n+1-[1:n]')
%    Local minima (several)
%    Search domain: -n^2 <= x <= n^2   (Montaz et al. (2005))
%    Cases considered: n = 10 Montaz et al. (2005)
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
dim  = size(x,1);
xaux = x(1:dim-1,1);
f    = sum((x-1).^2,1)-sum(x(2:dim,1).*xaux,1);
%
% End of neumaier3.

