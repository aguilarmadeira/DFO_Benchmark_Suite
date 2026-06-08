function varargout = fifteenn_local_minima_10D(varargin)
%FIFTEENN_LOCAL_MINIMA_10D  fifteenn_local_minima 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-62441.2107599, 18732.363228]   (range: 81173.5739879)
%   x2   ∈ [-37376.8092004, 11213.0427601]   (range: 48589.8519606)
%   x3   ∈ [-77482.9231597, 23244.8769479]   (range: 100727.800108)
%   x4   ∈ [-2410.82742815, 723.248228445]   (range: 3134.0756566)
%   x5   ∈ [-66301.479357, 19890.4438071]   (range: 86191.9231641)
%   x6   ∈ [-413221578.761, 123966473.628]   (range: 537188052.389)
%   x7   ∈ [-93835.4574004, 28150.6372201]   (range: 121986.09462)
%   x8   ∈ [-188736162.782, 56620848.8347]   (range: 245357011.617)
%   x9   ∈ [-55703.9358696, 16711.1807609]   (range: 72415.1166305)
%   x10  ∈ [-308343545.868, 92503063.7604]   (range: 400846609.629)
%
% Effective contrast ratio (max range / min range): 171402.388216
%
% Known global minimum (WORK-space):
%   x* = [6244.121075989628;3737.680920043807;7748.292315967134;241.0827428151256;6630.147935703601;41322157.87606096;9383.545740035828;18873616.27823228;5570.393586960949;30834354.58681077]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = fifteenn_local_minima_10D(n)   % Get bounds for dimension n (must be 10)
%   info = fifteenn_local_minima_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3;3;3;3;3];
lb_work = [-62441.21075989627;-37376.80920043807;-77482.92315967128;-2410.827428151253;-66301.47935703593;-413221578.7606097;-93835.45740035821;-188736162.7823228;-55703.93586960938;-308343545.8681074];
ub_work = [18732.36322796888;11213.04276013142;23244.87694790138;723.2482284453761;19890.44380711078;123966473.6281829;28150.63722010747;56620848.83469684;16711.18076088282;92503063.76043221];
scale_factors = [6244.121075989627;3737.680920043807;7748.292315967127;241.0827428151253;6630.147935703593;41322157.87606097;9383.545740035821;18873616.27823228;5570.393586960939;30834354.58681074];
contrast_ratio = 171402.388215522;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 19;
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
    info.x_global_min_orig = [1;1;1;1;1;1;1;1;1;1];
    info.x_global_min_work = [6244.121075989628;3737.680920043807;7748.292315967134;241.0827428151256;6630.147935703601;41322157.87606096;9383.545740035828;18873616.27823228;5570.393586960949;30834354.58681077];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=10): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = fifteenn_local_minima_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = fifteenn_local_minima_orig(x);
%
% Purpose:
%
%    Function fifteenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (15^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
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
f   = 0.1*sum(((x(1:dim-1)-1).^2).*(1+10*sin(3*pi.*x(2:dim)).^2),1);
f   = f + 0.1*(sin(3*pi*x(1))^2+(x(dim)-1)^2*(1+sin(2*pi*x(dim))^2));
%
% End of fifteenn_local_minima.

