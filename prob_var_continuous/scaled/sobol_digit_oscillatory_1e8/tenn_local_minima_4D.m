function varargout = tenn_local_minima_4D(varargin)
%TENN_LOCAL_MINIMA_4D  tenn_local_minima 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-569332930.256, 170799879.077]   (range: 740132809.333)
%   x2   ∈ [-366651785.072, 109995535.522]   (range: 476647320.594)
%   x3   ∈ [-941341072.006, 282402321.602]   (range: 1223743393.61)
%   x4   ∈ [-236879343.914, 71063803.1743]   (range: 307943147.088)
%
% Effective contrast ratio (max range / min range): 3.97392637303
%
% Known global minimum (WORK-space):
%   x* = [56933293.02559173;36665178.5072394;94134107.20057964;23687934.39141896]
%   f* = 0
%
% USAGE:
%   f = tenn_local_minima_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = tenn_local_minima_4D(n)   % Get bounds for dimension n (must be 4)
%   info = tenn_local_minima_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [3;3;3;3];
lb_work = [-569332930.255918;-366651785.0723941;-941341072.0057952;-236879343.9141897];
ub_work = [170799879.0767754;109995535.5217182;282402321.6017386;71063803.17425691];
scale_factors = [56933293.0255918;36665178.50723941;94134107.20057952;23687934.39141897];
contrast_ratio = 3.973926373026426;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'tenn_local_minima';
    info.source_p = 58;
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
    info.x_global_min_orig = [1;1;1;1];
    info.x_global_min_work = [56933293.02559173;36665178.5072394;94134107.20057964;23687934.39141896];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Ten Local Minima (n=4): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 4D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 4 components.');
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

