function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   schwefel (source instance p=46)
% Dimension: n = 10
% Strategy folder: sobol_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 500
% Effective contrast: 680.5457476716454
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = schwefel_orig(x_orig)

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [500;500;500;500;500;500;500;500;500;500];
lb_work = [-74644531.25;-324394.53125;-199519531.25;-449269.53125;-43425781.25;-293175.78125;-168300781.25;-418050.78125;-105863281.25;-355613.28125];
ub_work = [74644531.25;324394.53125;199519531.25;449269.53125;43425781.25;293175.78125;168300781.25;418050.78125;105863281.25;355613.28125];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

if nargin == 0
    info.name = mfilename;
    info.problem = 'schwefel';
    info.source_p = 46;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 500;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -4189.829;
    info.x_global_min_orig = [420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687];
    info.x_global_min_work = [62846022.56484374;273119.8882148437;167982955.3898438;378256.8210398437;36561789.35859375;246835.6550085937;141698722.1835938;351972.5878335937;89130255.77109376;299404.1214210937];
    info.global_min_note = 'Schwefel (n=10): x*=420.9687, f*=-418.9829*n. Ref: Ali et al. (2005).';
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
varargout{1} = schwefel_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schwefel_orig(x);
%
% Purpose:
%
%    Function schwefel is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -418.9829*n
%    Global minima = 420.97*ones(n,1)
%    Local minima 
%    Search domain: -500 <= x <= 500 (Montaz et al. (2005))
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
f = -sum(x.*sin(sqrt(abs(x))),1);
%
% End of schwefel.

