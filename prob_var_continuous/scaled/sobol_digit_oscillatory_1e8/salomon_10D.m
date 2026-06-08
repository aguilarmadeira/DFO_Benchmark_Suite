function varargout = salomon_10D(varargin)
%SALOMON_10D  salomon 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-5245865871.49, 1573759761.45]   (range: 6819625632.93)
%   x2   ∈ [-2383947599.91, 715184279.972]   (range: 3099131879.88)
%   x3   ∈ [-753507.070048, 226052.121014]   (range: 979559.191062)
%   x4   ∈ [-483182.856319, 144954.856896]   (range: 628137.713215)
%   x5   ∈ [-680923.64199, 204277.092597]   (range: 885200.734587)
%   x6   ∈ [-82194.1522339, 24658.2456702]   (range: 106852.397904)
%   x7   ∈ [-9097209948.83, 2729162984.65]   (range: 11826372933.5)
%   x8   ∈ [-3268591185.5, 980577355.651]   (range: 4249168541.16)
%   x9   ∈ [-600724.533384, 180217.360015]   (range: 780941.893399)
%   x10  ∈ [-135038.447621, 40511.5342863]   (range: 175549.981907)
%
% Effective contrast ratio (max range / min range): 110679.527699
%
% Known global minimum (WORK-space):
%   x* = [0;4.76837158203125e-07;0;0;1.164153218269348e-10;0;0;4.76837158203125e-07;1.164153218269348e-10;0]
%   f* = 0
%
% USAGE:
%   f = salomon_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = salomon_10D(n)   % Get bounds for dimension n (must be 10)
%   info = salomon_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [30;30;30;30;30;30;30;30;30;30];
lb_work = [-5245865871.486582;-2383947599.906651;-753507.0700476726;-483182.8563191153;-680923.6419897913;-82194.15223386715;-9097209948.832882;-3268591185.50434;-600724.5333837565;-135038.4476211536];
ub_work = [1573759761.445974;715184279.9719952;226052.1210143018;144954.8568957346;204277.0925969374;24658.24567016015;2729162984.649865;980577355.6513021;180217.360015127;40511.53428634606];
scale_factors = [52458658.71486582;23839475.99906651;7535.070700476726;4831.828563191152;6809.236419897913;821.9415223386715;90972099.48832883;32685911.8550434;6007.245333837565;1350.384476211535];
contrast_ratio = 110679.5276986199;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'salomon';
    info.source_p = 43;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;4.76837158203125e-07;0;0;1.164153218269348e-10;0;0;4.76837158203125e-07;1.164153218269348e-10;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Salomon (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = salomon_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x);
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
%
% End of salomon.

