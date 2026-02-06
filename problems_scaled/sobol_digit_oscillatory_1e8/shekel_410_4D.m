function varargout = shekel_410_4D(varargin)
%SHEKEL_410_4D  shekel_410 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 71332.2105894]   (range: 71332.2105894)
%   x2   ∈ [0           , 8395.44211728]   (range: 8395.44211728)
%   x3   ∈ [0           , 91189.5226637]   (range: 91189.5226637)
%   x4   ∈ [0           , 291377200.392]   (range: 291377200.392)
%
% Effective contrast ratio (max range / min range): 34706.5939258
%
% Known global minimum (WORK-space):
%   x* = [28532.88423574718;3358.176846913319;36475.80906548744;116550880.1569774]
%   f* = -10.5364
%
% USAGE:
%   f = shekel_410_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = shekel_410_4D(n)   % Get bounds for dimension n (must be 4)
%   info = shekel_410_4D()        % Get complete problem information

nloc = 4;
lb_orig = [0;0;0;0];
ub_orig = [10;10;10;10];
lb_work = [0;0;0;0];
ub_work = [71332.21058936795;8395.442117283295;91189.5226637186;291377200.3924434];
scale_factors = [7133.221058936795;839.5442117283296;9118.95226637186;29137720.03924434];
contrast_ratio = 34706.59392583972;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'shekel_410';
    info.source_p = 49;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -10.5364;
    info.x_global_min_orig = [4;4;4;4];
    info.x_global_min_work = [28532.88423574718;3358.176846913319;36475.80906548744;116550880.1569774];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Shekel-4,10 (4D): x*=(4,4,4,4), f*=-10.5364. Ref: Brachetti et al. (1997).';
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
varargout{1} = shekel_410_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = shekel_410_orig(x);
%
% Purpose:
%
%    Function shekel_4 is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value  
%    Global minima = 4*ones(4,1)
%    Local minima (m local minima at x = A(:,1:4))
%    Search domain: 0 <= x <= 10 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
%    Cases considered: m = 5 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 7 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 10 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
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
m = 10;
A = [4 4 4 4 0.1
    1 1 1 1 0.2
    8 8 8 8 0.2
    6 6 6 6 0.4
    3 7 3 7 0.4
    2 9 2 9 0.6
    5 5 3 3 0.3
    8 1 8 1 0.7
    6 2 6 2 0.5
    7 3.6 7 3.6 0.5];
f = -sum(1./(sum((repmat(x',m,1)-A(1:m,1:4)).^2,2)+A(1:m,5)),1);
%
% End of shekel_4.

