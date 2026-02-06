function varargout = kowalik_4D(varargin)
%KOWALIK_4D  kowalik 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 31.0932932021]   (range: 31.0932932021)
%   x2   ∈ [0           , 40515389.5821]   (range: 40515389.5821)
%   x3   ∈ [0           , 1802.21031673]   (range: 1802.21031673)
%   x4   ∈ [0           , 2542.62284242]   (range: 2542.62284242)
%
% Effective contrast ratio (max range / min range): 1303026.6469
%
% Known global minimum (WORK-space):
%   x* = [14.21407689238524;18328390.52525159;527.7901641857209;817.2716279204159]
%   f* = 0.00030748
%
% USAGE:
%   f = kowalik_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = kowalik_4D(n)   % Get bounds for dimension n (must be 4)
%   info = kowalik_4D()        % Get complete problem information

nloc = 4;
lb_orig = [0;0;0;0];
ub_orig = [0.42;0.42;0.42;0.42];
lb_work = [0;0;0;0];
ub_work = [31.09329320209271;40515389.58213509;1802.21031673173;2542.622842419071];
scale_factors = [74.03165048117313;96465213.29079783;4290.976944599357;6053.863910521598];
contrast_ratio = 1303026.646897866;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'kowalik';
    info.source_p = 28;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0.00030748;
    info.x_global_min_orig = [0.192;0.19;0.123;0.135];
    info.x_global_min_work = [14.21407689238524;18328390.52525159;527.7901641857209;817.2716279204159];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Kowalik (4D): f*=3.0748e-4. Ref: Ali et al. (2005).';
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
varargout{1} = kowalik_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = kowalik_orig(x);
%
% Purpose:
%
%    Function kowalik is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 3.0748*10^-4
%    Global minima = [0.192 0.190 0.123 0.135]'
%    Local minima 
%    Search domain: 0 <= x <= 0.42 (Montaz et al. (2005))
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
a = [0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b = [0.25 0.5 1 2 4 6 8 10 12 14 16];
f = sum((a-x(1).*(1+x(2).*b)./(1+x(3).*b+x(4).*b.^2)).^2,2);
%
% End of kowalik.

