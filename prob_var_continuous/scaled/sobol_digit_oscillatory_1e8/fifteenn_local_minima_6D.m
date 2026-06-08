function varargout = fifteenn_local_minima_6D(varargin)
%FIFTEENN_LOCAL_MINIMA_6D  fifteenn_local_minima 6D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-679145314.485, 203743594.346]   (range: 882888908.831)
%   x2   ∈ [-40050.3245179, 12015.0973554]   (range: 52065.4218732)
%   x3   ∈ [-80802.5022868, 24240.750686]   (range: 105043.252973)
%   x4   ∈ [-2158.74074065, 647.622222195]   (range: 2806.36296285)
%   x5   ∈ [-53398.7894197, 16019.6368259]   (range: 69418.4262457)
%   x6   ∈ [-26535.2199917, 7960.56599751]   (range: 34495.7859892)
%
% Effective contrast ratio (max range / min range): 314602.537348
%
% Known global minimum (WORK-space):
%   x* = [67914531.44854069;4005.032451785926;8080.250228675359;215.8740740651069;5339.878941973147;2653.521999169992]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_6D(x)          % Evaluate function at point x (6D vector)
%   [lb, ub] = fifteenn_local_minima_6D(n)   % Get bounds for dimension n (must be 6)
%   info = fifteenn_local_minima_6D()        % Get complete problem information

nloc = 6;
lb_orig = [-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3];
lb_work = [-679145314.4854078;-40050.32451785928;-80802.50228675368;-2158.740740651067;-53398.78941973147;-26535.21999169991];
ub_work = [203743594.3456224;12015.09735535779;24240.7506860261;647.6222221953202;16019.63682591944;7960.565997509973];
scale_factors = [67914531.44854078;4005.032451785929;8080.250228675367;215.8740740651068;5339.878941973147;2653.521999169991];
contrast_ratio = 314602.5373480376;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 17;
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
    info.x_global_min_work = [67914531.44854069;4005.032451785926;8080.250228675359;215.8740740651069;5339.878941973147;2653.521999169992];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=6): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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

