function varargout = rastrigin_10D(varargin)
%RASTRIGIN_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   rastrigin (source instance p=40)
% Dimension: n = 10
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 5.12
% Effective contrast: 68594.81975295641
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = rastrigin_orig(x_orig)

nloc = 10;
lb_orig = [-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12];
ub_orig = [5.12;5.12;5.12;5.12;5.12;5.12;5.12;5.12;5.12;5.12];
lb_work = [-38408.83789489685;-242634720.9046384;-379969000.6754542;-9536108.87159312;-459607952.0075021;-185105041.2217439;-30451.19299322746;-6700.330340727998;-418399986.8178775;-20895.53292350427];
ub_work = [38408.83789489685;242634720.9046384;379969000.6754542;9536108.87159312;459607952.0075021;185105041.2217439;30451.19299322746;6700.330340727998;418399986.8178775;20895.53292350427];
scale_factors = [7501.726151347041;47389593.9266872;74212695.44442464;1862521.263983031;89767178.12646525;36153328.36362185;5947.498631489738;1308.658269673437;81718747.4253667;4081.158774121928];
contrast_ratio = 68594.81975295641;

if nargin == 0
    info.name = mfilename;
    info.problem = 'rastrigin';
    info.source_p = 40;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Rastrigin (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = rastrigin_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = rastrigin_orig(x);
%
% Purpose:
%
%    Function rastrigin is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -5.12 <= x <= 5.12 (Montaz et al. (2005))
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
dim = size(x,1);
f   = 10*dim + sum(x.^2-10.*cos(2*pi.*x),1);
%
% End of rastrigin.

