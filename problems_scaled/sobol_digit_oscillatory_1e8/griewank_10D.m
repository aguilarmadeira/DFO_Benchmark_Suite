function varargout = griewank_10D(varargin)
%GRIEWANK_10D  griewank 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-2878157.7876, 2878157.7876]   (range: 5756315.5752)
%   x2   ∈ [-1866638.50396, 1866638.50396]   (range: 3733277.00792)
%   x3   ∈ [-3804917.8501, 3804917.8501]   (range: 7609835.7002)
%   x4   ∈ [-9495886763.09, 9495886763.09]   (range: 18991773526.2)
%   x5   ∈ [-21173907038.5, 21173907038.5]   (range: 42347814077 )
%   x6   ∈ [-1137148.78366, 1137148.78366]   (range: 2274297.56732)
%   x7   ∈ [-3190673.92312, 3190673.92312]   (range: 6381347.84625)
%   x8   ∈ [-541534223.819, 541534223.819]   (range: 1083068447.64)
%   x9   ∈ [-26094983527.1, 26094983527.1]   (range: 52189967054.2)
%   x10  ∈ [-16524474697.7, 16524474697.7]   (range: 33048949395.4)
%
% Effective contrast ratio (max range / min range): 22947.7302373
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = griewank_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = griewank_10D(n)   % Get bounds for dimension n (must be 10)
%   info = griewank_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-400;-400;-400;-400;-400;-400;-400;-400;-400;-400];
ub_orig = [400;400;400;400;400;400;400;400;400;400];
lb_work = [-2878157.787601171;-1866638.503958991;-3804917.850099366;-9495886763.086765;-21173907038.5177;-1137148.783662442;-3190673.923122885;-541534223.8191428;-26094983527.11537;-16524474697.67591];
ub_work = [2878157.787601171;1866638.503958991;3804917.850099366;9495886763.086765;21173907038.5177;1137148.783662442;3190673.923122885;541534223.8191428;26094983527.11537;16524474697.67591];
scale_factors = [7195.394469002928;4666.596259897477;9512.294625248414;23739716.90771692;52934767.59629425;2842.871959156105;7976.684807807213;1353835.559547857;65237458.81778843;41311186.74418976];
contrast_ratio = 22947.73023726116;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 22;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 400;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Griewank (n=10): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
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

