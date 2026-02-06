function varargout = neumaier3_10D(varargin)
%NEUMAIER3_10D  neumaier3 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-8569539805.66, 8569539805.66]   (range: 17139079611.3)
%   x2   ∈ [-336047.253108, 336047.253108]   (range: 672094.506216)
%   x3   ∈ [-5432900607.88, 5432900607.88]   (range: 10865801215.8)
%   x4   ∈ [-25253.5761669, 25253.5761669]   (range: 50507.1523337)
%   x5   ∈ [-9861567688.26, 9861567688.26]   (range: 19723135376.5)
%   x6   ∈ [-441799.478829, 441799.478829]   (range: 883598.957659)
%   x7   ∈ [-6720124389.38, 6720124389.38]   (range: 13440248778.8)
%   x8   ∈ [-130525.451378, 130525.451378]   (range: 261050.902755)
%   x9   ∈ [-778982.810635, 778982.810635]   (range: 1557965.62127)
%   x10  ∈ [-289305.789971, 289305.789971]   (range: 578611.579943)
%
% Effective contrast ratio (max range / min range): 390501.82925
%
% Known global minimum (WORK-space):
%   x* = [856953980.5662489;60488.50555947487;1303896145.89078;7071.001326719685;2958470306.478094;132539.8436488358;1881634829.027754;31326.1083306353;140216.905914303;28930.57899713988]
%   f* = -210
%
% USAGE:
%   f = neumaier3_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = neumaier3_10D(n)   % Get bounds for dimension n (must be 10)
%   info = neumaier3_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [100;100;100;100;100;100;100;100;100;100];
lb_work = [-8569539805.662486;-336047.2531081937;-5432900607.878252;-25253.57616685602;-9861567688.260309;-441799.4788294528;-6720124389.384836;-130525.4513776471;-778982.8106350169;-289305.7899713986];
ub_work = [8569539805.662486;336047.2531081937;5432900607.878252;25253.57616685602;9861567688.260309;441799.4788294528;6720124389.384836;130525.4513776471;778982.8106350169;289305.7899713986];
scale_factors = [85695398.05662486;3360.472531081938;54329006.07878252;252.5357616685602;98615676.88260309;4417.994788294528;67201243.89384836;1305.254513776471;7789.828106350169;2893.057899713986];
contrast_ratio = 390501.8292499537;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

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
    info.x_global_min_work = [856953980.5662489;60488.50555947487;1303896145.89078;7071.001326719685;2958470306.478094;132539.8436488358;1881634829.027754;31326.1083306353;140216.905914303;28930.57899713988];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Neumaier3 (10D): x*_i=i*(n+1-i), f*=-210. Ref: Ali et al. (2005).';
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

