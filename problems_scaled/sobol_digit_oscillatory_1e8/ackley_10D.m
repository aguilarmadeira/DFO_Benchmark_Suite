function varargout = ackley_10D(varargin)
%ACKLEY_10D  ackley 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-1757068234.52, 1757068234.52]   (range: 3514136469.04)
%   x2   ∈ [-105438.374236, 105438.374236]   (range: 210876.748472)
%   x3   ∈ [-2449343645.81, 2449343645.81]   (range: 4898687291.61)
%   x4   ∈ [-34128.2301534, 34128.2301534]   (range: 68256.4603068)
%   x5   ∈ [-200490.649095, 200490.649095]   (range: 400981.298191)
%   x6   ∈ [-120774.670881, 120774.670881]   (range: 241549.341762)
%   x7   ∈ [-281248.117387, 281248.117387]   (range: 562496.234774)
%   x8   ∈ [-420998235.451, 420998235.451]   (range: 841996470.902)
%   x9   ∈ [-154096.12599, 154096.12599]   (range: 308192.25198)
%   x10  ∈ [-93054.3748343, 93054.3748343]   (range: 186108.749669)
%
% Effective contrast ratio (max range / min range): 71768.844584
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
ub_orig = [30;30;30;30;30;30;30;30;30;30];
lb_work = [-1757068234.521907;-105438.3742357855;-2449343645.807012;-34128.23015341543;-200490.6490954587;-120774.6708810126;-281248.1173867774;-420998235.4510967;-154096.1259901716;-93054.37483433889];
ub_work = [1757068234.521907;105438.3742357855;2449343645.807012;34128.23015341543;200490.6490954587;120774.6708810126;281248.1173867774;420998235.4510967;154096.1259901716;93054.37483433889];
scale_factors = [58568941.15073022;3514.612474526182;81644788.19356707;1137.607671780514;6683.021636515289;4025.822362700419;9374.937246225913;14033274.51503656;5136.537533005721;3101.812494477963];
contrast_ratio = 71768.84458398703;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

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

