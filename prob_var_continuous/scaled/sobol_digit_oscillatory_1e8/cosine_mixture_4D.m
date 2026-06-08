function varargout = cosine_mixture_4D(varargin)
%COSINE_MIXTURE_4D  cosine_mixture 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-5555.18211359, 1666.55463408]   (range: 7221.73674767)
%   x2   ∈ [-18732753.5438, 5619826.06313]   (range: 24352579.6069)
%   x3   ∈ [-96576150.6673, 28972845.2002]   (range: 125548995.867)
%   x4   ∈ [-3395.31328139, 1018.59398442]   (range: 4413.90726581)
%
% Effective contrast ratio (max range / min range): 28443.9586758
%
% USAGE:
%   f = cosine_mixture_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = cosine_mixture_4D(n)   % Get bounds for dimension n (must be 4)
%   info = cosine_mixture_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-1;-1;-1;-1];
ub_orig = [0.3;0.3;0.3;0.3];
lb_work = [-5555.182113590938;-18732753.54378139;-96576150.66727184;-3395.313281393187];
ub_work = [1666.554634077282;5619826.063134417;28972845.20018156;1018.593984417956];
scale_factors = [5555.182113590938;18732753.54378139;96576150.66727184;3395.313281393187];
contrast_ratio = 28443.95867577912;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'cosine_mixture';
    info.source_p = 9;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 1;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = false;
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
varargout{1} = cosine_mixture_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = cosine_mixture_orig(x);
%
% Purpose:
%
%    Function cosine_mixture is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value: -2 if n = 2
%                          -4 if n = 4
%    Global minima
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
f = 0.1*sum(cos(5*pi.*x),1) - sum(x.^2,1); 
%
% End of cosine_mixture.

