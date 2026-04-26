function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  schwefel 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-4125100.6554, 1237530.19662]   (range: 5362630.85203)
%   x2   ∈ [-344709.323531, 103412.797059]   (range: 448122.120591)
%   x3   ∈ [-3109097.09437, 932729.128312]   (range: 4041826.22269)
%   x4   ∈ [-18312455772.8, 5493736731.83]   (range: 23806192504.6)
%   x5   ∈ [-49162911497.1, 14748873449.1]   (range: 63911784946.2)
%   x6   ∈ [-1116568.37666, 334970.512998]   (range: 1451538.88966)
%   x7   ∈ [-3587251.60849, 1076175.48255]   (range: 4663427.09103)
%   x8   ∈ [-22899828275.3, 6869948482.59]   (range: 29769776757.9)
%   x9   ∈ [-3866221.1065, 1159866.33195]   (range: 5026087.43846)
%   x10  ∈ [-144725.788217, 43417.7364652]   (range: 188143.524682)
%
% Effective contrast ratio (max range / min range): 339696.968333
%
% USAGE:
%   f = schwefel_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = schwefel_10D(n)   % Get bounds for dimension n (must be 10)
%   info = schwefel_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [150;150;150;150;150;150;150;150;150;150];
lb_work = [-4125100.655404098;-344709.323531174;-3109097.094374723;-18312455772.78201;-49162911497.05669;-1116568.376659427;-3587251.60848704;-22899828275.29753;-3866221.10650436;-144725.7882172262];
ub_work = [1237530.196621229;103412.7970593522;932729.1283124168;5493736731.834602;14748873449.117;334970.5129978281;1076175.482546112;6869948482.589258;1159866.331951308;43417.73646516787];
scale_factors = [8250.201310808196;689.418647062348;6218.194188749446;36624911.54556401;98325822.99411337;2233.136753318854;7174.50321697408;45799656.55059506;7732.442213008721;289.4515764344525];
contrast_ratio = 339696.9683334223;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

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

