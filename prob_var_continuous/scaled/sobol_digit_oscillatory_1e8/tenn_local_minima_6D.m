function varargout = tenn_local_minima_6D(varargin)
%TENN_LOCAL_MINIMA_6D  tenn_local_minima 6D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-815592352.663, 244677705.799]   (range: 1060270058.46)
%   x2   ∈ [-19558506.0456, 5867551.81368]   (range: 25426057.8593)
%   x3   ∈ [-658889457.482, 197666837.245]   (range: 856556294.726)
%   x4   ∈ [-48675.436544, 14602.6309632]   (range: 63278.0675072)
%   x5   ∈ [-93790.7108227, 28137.2132468]   (range: 121927.924069)
%   x6   ∈ [-141327569.531, 42398270.8592]   (range: 183725840.39)
%
% Effective contrast ratio (max range / min range): 16755.7275408
%
% Known global minimum (WORK-space):
%   x* = [81559235.2662729;1955850.604561552;65888945.74818492;4867.543654399728;9379.071082267707;14132756.95305681]
%   f* = 0
%
% USAGE:
%   f = tenn_local_minima_6D(x)          % Evaluate function at point x (6D vector)
%   [lb, ub] = tenn_local_minima_6D(n)   % Get bounds for dimension n (must be 6)
%   info = tenn_local_minima_6D()        % Get complete problem information

nloc = 6;
lb_orig = [-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3];
lb_work = [-815592352.6627285;-19558506.04561551;-658889457.4818492;-48675.43654399729;-93790.71082267709;-141327569.5305684];
ub_work = [244677705.7988186;5867551.813684653;197666837.2445548;14602.63096319919;28137.21324680313;42398270.85917051];
scale_factors = [81559235.26627286;1955850.604561551;65888945.74818492;4867.543654399728;9379.07108226771;14132756.95305683];
contrast_ratio = 16755.7275408413;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'tenn_local_minima';
    info.source_p = 59;
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
    info.x_global_min_orig = [1;1;1;1;1;1];
    info.x_global_min_work = [81559235.2662729;1955850.604561552;65888945.74818492;4867.543654399728;9379.071082267707;14132756.95305681];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Ten Local Minima (n=6): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 6D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 6 components.');
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

