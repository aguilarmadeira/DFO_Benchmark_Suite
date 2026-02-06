function varargout = goldstein_price_2D(varargin)
%GOLDSTEIN_PRICE_2D  goldstein_price 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (BASELINE HETEROGENEITY):
%
%   x1   ∈ [-2          , 2           ]   (range: 4           )
%   x2   ∈ [-2          , 2           ]   (range: 4           )
%
% Effective contrast ratio (max range / min range): 1
%
% Known global minimum (WORK-space):
%   x* = [0;-1]
%   f* = 3
%
% USAGE:
%   f = goldstein_price_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = goldstein_price_2D(n)   % Get bounds for dimension n (must be 2)
%   info = goldstein_price_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-2;-2];
ub_orig = [2;2];
lb_work = [-2;-2];
ub_work = [2;2];
scale_factors = [1;1];
contrast_ratio = 1;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'goldstein_price';
    info.source_p = 21;
    info.dimension = nloc;
    info.strategy = 'baseline';
    info.kappa = 1;
    info.bound_p = 2;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 3;
    info.x_global_min_orig = [0;-1];
    info.x_global_min_work = [0;-1];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Goldstein-Price (2D): x*=(0,-1), f*=3. Ref: Brachetti et al. (1997).';
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
varargout{1} = goldstein_price_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = goldstein_price_orig(x);
%
% Purpose:
%
%    Function goldstein_price is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = 3
%    Global minima = [0 -1]'
%    Local minima (four)
%    Search domain: -2 <= x <= 2 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
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
f = 1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2);
f = f * (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2));
%
% End of goldstein_price.

