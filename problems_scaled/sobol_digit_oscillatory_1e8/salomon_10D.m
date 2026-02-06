function varargout = salomon_10D(varargin)
%SALOMON_10D  salomon 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-246629.611888, 246629.611888]   (range: 493259.223777)
%   x2   ∈ [-9072935261.86, 9072935261.86]   (range: 18145870523.7)
%   x3   ∈ [-2776159306.8, 2776159306.8]   (range: 5552318613.59)
%   x4   ∈ [-5637788435.77, 5637788435.77]   (range: 11275576871.5)
%   x5   ∈ [-68133.4310978, 68133.4310978]   (range: 136266.862196)
%   x6   ∈ [-774299.090579, 774299.090579]   (range: 1548598.18116)
%   x7   ∈ [-4120481211.6, 4120481211.6]   (range: 8240962423.2)
%   x8   ∈ [-7427897223.84, 7427897223.84]   (range: 14855794447.7)
%   x9   ∈ [-173786.342564, 173786.342564]   (range: 347572.685127)
%   x10  ∈ [-9811144539.59, 9811144539.59]   (range: 19622289079.2)
%
% Effective contrast ratio (max range / min range): 143998.979378
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = salomon_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = salomon_10D(n)   % Get bounds for dimension n (must be 10)
%   info = salomon_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [100;100;100;100;100;100;100;100;100;100];
lb_work = [-246629.6118884623;-9072935261.855909;-2776159306.797195;-5637788435.765785;-68133.43109778507;-774299.0905786331;-4120481211.602268;-7427897223.844584;-173786.3425635711;-9811144539.58934];
ub_work = [246629.6118884623;9072935261.855909;2776159306.797195;5637788435.765785;68133.43109778507;774299.0905786331;4120481211.602268;7427897223.844584;173786.3425635711;9811144539.58934];
scale_factors = [2466.296118884623;90729352.61855909;27761593.06797195;56377884.35765785;681.3343109778507;7742.990905786331;41204812.11602268;74278972.23844583;1737.863425635712;98111445.39589339];
contrast_ratio = 143998.9793778092;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'salomon';
    info.source_p = 43;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Salomon (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = salomon_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x);
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
%
% End of salomon.

