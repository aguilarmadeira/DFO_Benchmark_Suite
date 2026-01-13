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
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 500
% Effective contrast: 114576.0960287952
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = schwefel_orig(x_orig)

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [500;500;500;500;500;500;500;500;500;500];
lb_work = [-3302635.461151428;-422667.5957036936;-48427603033.60638;-14744455755.24443;-3029274.888182676;-774144.3876961047;-3866269.790242362;-23728200186.24872;-3566398.521413769;-1979645979.555758];
ub_work = [3302635.461151428;422667.5957036936;48427603033.60638;14744455755.24443;3029274.888182676;774144.3876961047;3866269.790242362;23728200186.24872;3566398.521413769;1979645979.555758];
scale_factors = [6605.270922302855;845.3351914073872;96855206.06721276;29488911.51048886;6058.549776365353;1548.288775392209;7732.539580484724;47456400.37249745;7132.797042827538;3959291.959111516];
contrast_ratio = 114576.0960287952;

if nargin == 0
    info.name = mfilename;
    info.problem = 'schwefel';
    info.source_p = 46;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 500;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -4189.829;
    info.x_global_min_orig = [420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687];
    info.x_global_min_work = [2780612.313309634;355859.656591019;40773010186.34666;12413908742.98553;2550459.823241813;651781.1130014503;3255157.134895199;19977659171.48977;3002684.298482953;1666737988.947628];
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

