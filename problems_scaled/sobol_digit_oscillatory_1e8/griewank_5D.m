function varargout = griewank_5D(varargin)
%GRIEWANK_5D  griewank 5D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-32552574667.1, 32552574667.1]   (range: 65105149334.2)
%   x2   ∈ [-26056456536, 26056456536 ]   (range: 52112913072.1)
%   x3   ∈ [-5973929.11153, 5973929.11153]   (range: 11947858.2231)
%   x4   ∈ [-6372659848.66, 6372659848.66]   (range: 12745319697.3)
%   x5   ∈ [-4282187.27975, 4282187.27975]   (range: 8564374.55951)
%
% Effective contrast ratio (max range / min range): 7601.85684102
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = griewank_5D(x)          % Evaluate function at point x (5D vector)
%   [lb, ub] = griewank_5D(n)   % Get bounds for dimension n (must be 5)
%   info = griewank_5D()        % Get complete problem information

nloc = 5;
lb_orig = [-600;-600;-600;-600;-600];
ub_orig = [600;600;600;600;600];
lb_work = [-32552574667.12121;-26056456536.02848;-5973929.111532643;-6372659848.658818;-4282187.279753974];
ub_work = [32552574667.12121;26056456536.02848;5973929.111532643;6372659848.658818;4282187.279753974];
scale_factors = [54254291.11186868;43427427.56004746;9956.548519221073;10621099.7477647;7136.978799589956];
contrast_ratio = 7601.856841018748;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 23;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 600;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Griewank (n=5): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 5D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 5 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = griewank_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = griewank_orig(x);
%
% Purpose:
%
%    Function griewank is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -600 <= x <= 600 (Huyer and Neumaier (1999))
%                   -400 <= x <= 400 (Storn and Price (1997))
%    Cases considered: n = 10 Storn and Price (1997)
%                      n = 5 Huyer and Neumaier (1999), ICEO
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
f = 1 + 1/4000*sum(x.^2,1)- prod(cos(x./sqrt([1:n]')),1);
%
% End of griewank.

