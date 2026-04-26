function varargout = neumaier3_10D(varargin)
%NEUMAIER3_10D  neumaier3 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-172328.507033, 51698.5521098]   (range: 224027.059142)
%   x2   ∈ [-8795680354.96, 2638704106.49]   (range: 11434384461.4)
%   x3   ∈ [-287782.6602, 86334.7980601]   (range: 374117.45826)
%   x4   ∈ [-517957.010576, 155387.103173]   (range: 673344.113749)
%   x5   ∈ [-20538.8123812, 6161.64371435]   (range: 26700.4560955)
%   x6   ∈ [-7814132232.55, 2344239669.76]   (range: 10158371902.3)
%   x7   ∈ [-377086.212051, 113125.863615]   (range: 490212.075667)
%   x8   ∈ [-6785979227.83, 2035793768.35]   (range: 8821772996.17)
%   x9   ∈ [-2346967747.74, 704090324.322]   (range: 3051058072.06)
%   x10  ∈ [-9419430975.31, 2825829292.59]   (range: 12245260267.9)
%
% Effective contrast ratio (max range / min range): 458616.145885
%
% Known global minimum (WORK-space):
%   x* = [17232.85070325813;1583222463.892746;69067.83844807482;145027.9629612433;6161.643714352107;2344239669.764274;105584.1393743902;1628635014.678417;422454194.5930567;941943097.531477]
%   f* = -210
%
% USAGE:
%   f = neumaier3_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = neumaier3_10D(n)   % Get bounds for dimension n (must be 10)
%   info = neumaier3_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [30;30;30;30;30;30;30;30;30;30];
lb_work = [-172328.5070325813;-8795680354.959703;-287782.6602003119;-517957.0105758691;-20538.81238117368;-7814132232.54758;-377086.2120513937;-6785979227.826735;-2346967747.739206;-9419430975.31477];
ub_work = [51698.55210977438;2638704106.487911;86334.79806009357;155387.1031727607;6161.643714352106;2344239669.764274;113125.8636154181;2035793768.34802;704090324.3217616;2825829292.594431];
scale_factors = [1723.285070325813;87956803.54959704;2877.826602003119;5179.570105758691;205.3881238117368;78141322.3254758;3770.862120513937;67859792.27826735;23469677.47739206;94194309.75314769];
contrast_ratio = 458616.1458852811;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

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
    info.x_global_min_work = [17232.85070325813;1583222463.892746;69067.83844807482;145027.9629612433;6161.643714352107;2344239669.764274;105584.1393743902;1628635014.678417;422454194.5930567;941943097.531477];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Some minimizer(s) lie on the ORIGINAL boundary (t≈0 or t≈1). With clip01 in decoding, WORK-space minimizers may be non-unique; returned x*_work is a canonical representative. Original note: Neumaier3 (10D): x*_i=i*(n+1-i), f*=-210. Ref: Ali et al. (2005).';
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

