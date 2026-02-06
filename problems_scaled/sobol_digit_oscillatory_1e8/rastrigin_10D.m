function varargout = rastrigin_10D(varargin)
%RASTRIGIN_10D  rastrigin 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-503143240.784, 503143240.784]   (range: 1006286481.57)
%   x2   ∈ [-233340688.068, 233340688.068]   (range: 466681376.137)
%   x3   ∈ [-269240403.556, 269240403.556]   (range: 538480807.111)
%   x4   ∈ [-19311823.5938, 19311823.5938]   (range: 38623647.1875)
%   x5   ∈ [-41466.442455, 41466.442455]   (range: 82932.8849101)
%   x6   ∈ [-128866385.865, 128866385.865]   (range: 257732771.731)
%   x7   ∈ [-356721513.135, 356721513.135]   (range: 713443026.269)
%   x8   ∈ [-122792324.788, 122792324.788]   (range: 245584649.577)
%   x9   ∈ [-46881.1772193, 46881.1772193]   (range: 93762.3544387)
%   x10  ∈ [-20276.0673932, 20276.0673932]   (range: 40552.1347863)
%
% Effective contrast ratio (max range / min range): 24814.6364395
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = rastrigin_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = rastrigin_10D(n)   % Get bounds for dimension n (must be 10)
%   info = rastrigin_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12];
ub_orig = [5.12;5.12;5.12;5.12;5.12;5.12;5.12;5.12;5.12;5.12];
lb_work = [-503143240.7844462;-233340688.0684597;-269240403.5557285;-19311823.59376011;-41466.4424550487;-128866385.8653758;-356721513.1346471;-122792324.7883925;-46881.17721933345;-20276.06739316135];
ub_work = [503143240.7844462;233340688.0684597;269240403.5557285;19311823.59376011;41466.4424550487;128866385.8653758;356721513.1346471;122792324.7883925;46881.17721933345;20276.06739316135];
scale_factors = [98270164.21571214;45574353.13837104;52586016.31947822;3771840.545656271;8098.914542001699;25169215.9893312;69672170.53411075;23982875.93523291;9156.479925651063;3960.169412726826];
contrast_ratio = 24814.6364395171;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

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
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Rastrigin (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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

