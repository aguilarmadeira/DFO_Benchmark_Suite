function varargout = ackley_10D(varargin)
%ACKLEY_10D  ackley 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-2528.9707677, 758.691230311]   (range: 3287.66199801)
%   x2   ∈ [-2577963610.77, 773389083.231]   (range: 3351352694  )
%   x3   ∈ [-1024015174.1, 307204552.231]   (range: 1331219726.33)
%   x4   ∈ [-167041.411229, 50112.4233688]   (range: 217153.834598)
%   x5   ∈ [-40551.524134, 12165.4572402]   (range: 52716.9813742)
%   x6   ∈ [-2945006877.52, 883502063.257]   (range: 3828508940.78)
%   x7   ∈ [-140534.067489, 42160.2202466]   (range: 182694.287735)
%   x8   ∈ [-204148.672393, 61244.601718]   (range: 265393.274112)
%   x9   ∈ [-233439280.955, 70031784.2866]   (range: 303471065.242)
%   x10  ∈ [-237017.877512, 71105.3632537]   (range: 308123.240766)
%
% Effective contrast ratio (max range / min range): 1164508.07385
%
% Known global minimum (WORK-space):
%   x* = [4.547473508864641e-13;0;0;0;0;0;0;2.91038304567337e-11;2.980232238769531e-08;2.91038304567337e-11]
%   f* = 0
%
% USAGE:
%   f = ackley_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = ackley_10D(n)   % Get bounds for dimension n (must be 10)
%   info = ackley_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-30;-30;-30;-30;-30;-30;-30;-30;-30;-30];
ub_orig = [9;9;9;9;9;9;9;9;9;9];
lb_work = [-2528.970767702465;-2577963610.770712;-1024015174.101778;-167041.4112292735;-40551.52413402626;-2945006877.522953;-140534.0674886398;-204148.6723934775;-233439280.9554268;-237017.877512333];
ub_work = [758.6912303107396;773389083.2312136;307204552.2305334;50112.42336878207;12165.45724020788;883502063.2568859;42160.22024659193;61244.60171804325;70031784.28662805;71105.36325369991];
scale_factors = [84.29902559008217;85932120.35902373;34133839.13672594;5568.047040975785;1351.717471134209;98166895.91743177;4684.468916287992;6804.955746449249;7781309.365180894;7900.595917077766];
contrast_ratio = 1164508.073851107;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'ackley';
    info.source_p = 1;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 30;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [4.547473508864641e-13;0;0;0;0;0;0;2.91038304567337e-11;2.980232238769531e-08;2.91038304567337e-11];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Ackley (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = ackley_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = ackley_orig(x);
%
% Purpose:
%
%    Function ackley is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -30 <= x <= 30 (Montaz et al. (2005))
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
n = size(x,1);
f = -20*exp(-0.02*sqrt(1/n*sum(x.^2,1)))-exp(1/n*sum(cos(2*pi.*x),1))+20+exp(1);
%
% End of ackley.

