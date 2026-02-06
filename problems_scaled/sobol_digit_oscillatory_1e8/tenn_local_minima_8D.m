function varargout = tenn_local_minima_8D(varargin)
%TENN_LOCAL_MINIMA_8D  tenn_local_minima 8D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-2146.65858316, 2146.65858316]   (range: 4293.31716633)
%   x2   ∈ [-68467.4399267, 68467.4399267]   (range: 136934.879853)
%   x3   ∈ [-33570.3855984, 33570.3855984]   (range: 67140.7711968)
%   x4   ∈ [-99928.2902175, 99928.2902175]   (range: 199856.580435)
%   x5   ∈ [-218631225.521, 218631225.521]   (range: 437262451.042)
%   x6   ∈ [-616488822.994, 616488822.994]   (range: 1232977645.99)
%   x7   ∈ [-40397.7334059, 40397.7334059]   (range: 80795.4668117)
%   x8   ∈ [-80219.7954567, 80219.7954567]   (range: 160439.590913)
%
% Effective contrast ratio (max range / min range): 287185.315741
%
% Known global minimum (WORK-space):
%   x* = [214.6658583163517;6846.743992673044;3357.038559840585;9992.829021746875;21863122.55211776;61648882.29939187;4039.773340585147;8021.979545669936]
%   f* = 0
%
% USAGE:
%   f = tenn_local_minima_8D(x)          % Evaluate function at point x (8D vector)
%   [lb, ub] = tenn_local_minima_8D(n)   % Get bounds for dimension n (must be 8)
%   info = tenn_local_minima_8D()        % Get complete problem information

nloc = 8;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10];
lb_work = [-2146.658583163516;-68467.43992673043;-33570.38559840578;-99928.2902174686;-218631225.5211775;-616488822.9939183;-40397.73340585146;-80219.79545669924];
ub_work = [2146.658583163516;68467.43992673043;33570.38559840578;99928.2902174686;218631225.5211775;616488822.9939183;40397.73340585146;80219.79545669924];
scale_factors = [214.6658583163516;6846.743992673042;3357.038559840578;9992.82902174686;21863122.55211775;61648882.29939183;4039.773340585146;8021.979545669925];
contrast_ratio = 287185.3157409889;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

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
    info.x_global_min_work = [214.6658583163517;6846.743992673044;3357.038559840585;9992.829021746875;21863122.55211776;61648882.29939187;4039.773340585147;8021.979545669936];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Ten Local Minima (n=8): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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

