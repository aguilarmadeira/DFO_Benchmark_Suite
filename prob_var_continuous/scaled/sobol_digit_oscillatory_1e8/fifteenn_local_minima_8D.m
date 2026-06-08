function varargout = fifteenn_local_minima_8D(varargin)
%FIFTEENN_LOCAL_MINIMA_8D  fifteenn_local_minima 8D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-549617431.337, 164885229.401]   (range: 714502660.738)
%   x2   ∈ [-22879.1510984, 6863.74532951]   (range: 29742.8964279)
%   x3   ∈ [-86644.5672507, 25993.3701752]   (range: 112637.937426)
%   x4   ∈ [-42173.8412618, 12652.1523785]   (range: 54825.9936403)
%   x5   ∈ [-62923.5024569, 18877.0507371]   (range: 81800.553194)
%   x6   ∈ [-90649823.8512, 27194947.1554]   (range: 117844771.007)
%   x7   ∈ [-93917.4619771, 28175.2385931]   (range: 122092.70057)
%   x8   ∈ [-274739875.732, 82421962.7195]   (range: 357161838.451)
%
% Effective contrast ratio (max range / min range): 24022.6321761
%
% Known global minimum (WORK-space):
%   x* = [54961743.13368249;2287.915109836809;8664.456725070544;4217.384126179313;6292.350245691443;9064982.385120869;9391.746197713743;27473987.5731585]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_8D(x)          % Evaluate function at point x (8D vector)
%   [lb, ub] = fifteenn_local_minima_8D(n)   % Get bounds for dimension n (must be 8)
%   info = fifteenn_local_minima_8D()        % Get complete problem information

nloc = 8;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3;3;3];
lb_work = [-549617431.3368258;-22879.15109836811;-86644.5672507055;-42173.84126179313;-62923.5024569144;-90649823.85120872;-93917.46197713746;-274739875.7315852];
ub_work = [164885229.4010478;6863.745329510432;25993.37017521165;12652.15237853794;18877.05073707433;27194947.15536261;28175.23859314124;82421962.71947557];
scale_factors = [54961743.13368258;2287.915109836811;8664.45672507055;4217.384126179313;6292.35024569144;9064982.385120872;9391.746197713745;27473987.57315852];
contrast_ratio = 24022.63217607003;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 18;
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
    info.x_global_min_work = [54961743.13368249;2287.915109836809;8664.456725070544;4217.384126179313;6292.350245691443;9064982.385120869;9391.746197713743;27473987.5731585];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=8): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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

