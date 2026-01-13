function varargout = tenn_local_minima_8D(varargin)
%TENN_LOCAL_MINIMA_8D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   tenn_local_minima (source instance p=60)
% Dimension: n = 8
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 10
% Effective contrast: 84650.6123038856
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = tenn_local_minima_orig(x_orig)

nloc = 8;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10];
lb_work = [-93174.22207509597;-185397701.1378443;-528752720.7038549;-276194028.3111725;-84183.63016281961;-8784.99790562978;-743655451.7999135;-498712219.7493586];
ub_work = [93174.22207509597;185397701.1378443;528752720.7038549;276194028.3111725;84183.63016281961;8784.99790562978;743655451.7999135;498712219.7493586];
scale_factors = [9317.422207509597;18539770.11378443;52875272.07038549;27619402.83111725;8418.363016281961;878.499790562978;74365545.17999135;49871221.97493586];
contrast_ratio = 84650.6123038856;

if nargin == 0
    info.name = mfilename;
    info.problem = 'tenn_local_minima';
    info.source_p = 60;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;1;1;1;1;1;1;1];
    info.x_global_min_work = [9317.422207509604;18539770.11378443;52875272.07038558;27619402.83111727;8418.363016281975;878.4997905629789;74365545.17999136;49871221.97493595];
    info.global_min_note = 'Ten Local Minima (n=8): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 8D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 8 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = tenn_local_minima_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = tenn_local_minima_orig(x);
%
% Purpose:
%
%    Function tenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (10^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
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
f   = pi/dim*sum(((x(1:dim-1)-1).^2).*(1+10*sin(pi.*x(2:dim)).^2),1);
f   = f + pi/dim*(10*sin(pi*x(1))^2+(x(dim)-1)^2);
%
% End of tenn_local_minima.

