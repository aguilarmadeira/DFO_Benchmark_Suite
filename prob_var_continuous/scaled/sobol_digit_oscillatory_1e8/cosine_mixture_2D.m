function varargout = cosine_mixture_2D(varargin)
%COSINE_MIXTURE_2D  cosine_mixture 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-69243327.7197, 20772998.3159]   (range: 90016326.0357)
%   x2   ∈ [-18890718.5375, 5667215.56125]   (range: 24557934.0988)
%
% Effective contrast ratio (max range / min range): 3.66546818123
%
% USAGE:
%   f = cosine_mixture_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = cosine_mixture_2D(n)   % Get bounds for dimension n (must be 2)
%   info = cosine_mixture_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-1;-1];
ub_orig = [0.3;0.3];
lb_work = [-69243327.71974449;-18890718.53751341];
ub_work = [20772998.31592335;5667215.561254022];
scale_factors = [69243327.71974449;18890718.53751341];
contrast_ratio = 3.665468181225627;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'cosine_mixture';
    info.source_p = 8;
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
    if arg1 ~= nloc, error('This instance is 2D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 2 components.');
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

