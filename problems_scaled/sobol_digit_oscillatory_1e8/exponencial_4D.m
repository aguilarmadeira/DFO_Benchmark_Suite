function varargout = exponencial_4D(varargin)
%EXPONENCIAL_4D  exponencial 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-90169985.7998, 90169985.7998]   (range: 180339971.6 )
%   x2   ∈ [-29101477.8791, 29101477.8791]   (range: 58202955.7582)
%   x3   ∈ [-70630521.8487, 70630521.8487]   (range: 141261043.697)
%   x4   ∈ [-9530687.21889, 9530687.21889]   (range: 19061374.4378)
%
% Effective contrast ratio (max range / min range): 9.46101616063
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0]
%   f* = -1
%
% USAGE:
%   f = exponencial_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = exponencial_4D(n)   % Get bounds for dimension n (must be 4)
%   info = exponencial_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-1;-1;-1;-1];
ub_orig = [1;1;1;1];
lb_work = [-90169985.79982643;-29101477.87908779;-70630521.84874322;-9530687.218893221];
ub_work = [90169985.79982643;29101477.87908779;70630521.84874322;9530687.218893221];
scale_factors = [90169985.79982643;29101477.87908779;70630521.84874322;9530687.218893221];
contrast_ratio = 9.461016160626629;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'exponencial';
    info.source_p = 14;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 1;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1;
    info.x_global_min_orig = [0;0;0;0];
    info.x_global_min_work = [0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Exponencial (n=4): x*=0, f*=-1. Ref: Brachetti et al. (1997).';
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
varargout{1} = exponencial_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = exponencial_orig(x);
%
% Purpose:
%
%    Function exponencial is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = -1
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -1 <= x <= 1 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
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
f = -exp(-0.5*sum(x.^2,1)); 
%
% End of exponencial.

