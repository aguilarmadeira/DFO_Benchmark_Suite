function varargout = neumaier3_10D(varargin)
%NEUMAIER3_10D  neumaier3 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-14928906.25, 4478671.875 ]   (range: 19407578.125)
%   x2   ∈ [-64878.90625, 19463.671875]   (range: 84342.578125)
%   x3   ∈ [-39903906.25, 11971171.875]   (range: 51875078.125)
%   x4   ∈ [-89853.90625, 26956.171875]   (range: 116810.078125)
%   x5   ∈ [-8685156.25 , 2605546.875 ]   (range: 11290703.125)
%   x6   ∈ [-58635.15625, 17590.546875]   (range: 76225.703125)
%   x7   ∈ [-33660156.25, 10098046.875]   (range: 43758203.125)
%   x8   ∈ [-83610.15625, 25083.046875]   (range: 108693.203125)
%   x9   ∈ [-21172656.25, 6351796.875 ]   (range: 27524453.125)
%   x10  ∈ [-71122.65625, 21336.796875]   (range: 92459.453125)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   x* = [1492890.625;11678.203125;9576937.5;25159.09375;2605546.875;17590.546875;9424843.75;20066.4375;3811078.125;7112.265625]
%   f* = -210
%
% USAGE:
%   f = neumaier3_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = neumaier3_10D(n)   % Get bounds for dimension n (must be 10)
%   info = neumaier3_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [30;30;30;30;30;30;30;30;30;30];
lb_work = [-14928906.25;-64878.90625;-39903906.25;-89853.90625;-8685156.25;-58635.15625;-33660156.25;-83610.15625;-21172656.25;-71122.65625];
ub_work = [4478671.875;19463.671875;11971171.875;26956.171875;2605546.875;17590.546875;10098046.875;25083.046875;6351796.875;21336.796875];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

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
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -210;
    info.x_global_min_orig = [10;18;24;28;30;30;28;24;18;10];
    info.x_global_min_work = [1492890.625;11678.203125;9576937.5;25159.09375;2605546.875;17590.546875;9424843.75;20066.4375;3811078.125;7112.265625];
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

