function varargout = ackley_10D(varargin)
%ACKLEY_10D  ackley 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-4478671.875, 1343601.5625]   (range: 5822273.4375)
%   x2   ∈ [-19463.671875, 5839.1015625]   (range: 25302.7734375)
%   x3   ∈ [-11971171.875, 3591351.5625]   (range: 15562523.4375)
%   x4   ∈ [-26956.171875, 8086.8515625]   (range: 35043.0234375)
%   x5   ∈ [-2605546.875, 781664.0625 ]   (range: 3387210.9375)
%   x6   ∈ [-17590.546875, 5277.1640625]   (range: 22867.7109375)
%   x7   ∈ [-10098046.875, 3029414.0625]   (range: 13127460.9375)
%   x8   ∈ [-25083.046875, 7524.9140625]   (range: 32607.9609375)
%   x9   ∈ [-6351796.875, 1905539.0625]   (range: 8257335.9375)
%   x10  ∈ [-21336.796875, 6401.0390625]   (range: 27737.8359375)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = ackley_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = ackley_10D(n)   % Get bounds for dimension n (must be 10)
%   info = ackley_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-30;-30;-30;-30;-30;-30;-30;-30;-30;-30];
ub_orig = [9;9;9;9;9;9;9;9;9;9];
lb_work = [-4478671.875;-19463.671875;-11971171.875;-26956.171875;-2605546.875;-17590.546875;-10098046.875;-25083.046875;-6351796.875;-21336.796875];
ub_work = [1343601.5625;5839.1015625;3591351.5625;8086.8515625;781664.0625;5277.1640625;3029414.0625;7524.9140625;1905539.0625;6401.0390625];
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
    info.problem = 'ackley';
    info.source_p = 1;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 30;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
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

