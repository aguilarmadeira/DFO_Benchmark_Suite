function varargout = transistor_9D(varargin)
%TRANSISTOR_9D  transistor 9D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-892907440.667, 267872232.2 ]   (range: 1160779672.87)
%   x2   ∈ [-325216597.072, 97564979.1217]   (range: 422781576.194)
%   x3   ∈ [-66650.0187577, 19995.0056273]   (range: 86645.024385)
%   x4   ∈ [-11418.6879406, 3425.60638218]   (range: 14844.2943228)
%   x5   ∈ [-837085610.24, 251125683.072]   (range: 1088211293.31)
%   x6   ∈ [-382419663.415, 114725899.024]   (range: 497145562.439)
%   x7   ∈ [-59481.8581998, 17844.55746 ]   (range: 77326.4156598)
%   x8   ∈ [-187265039.624, 56179511.8872]   (range: 243444551.511)
%   x9   ∈ [-981037281.223, 294311184.367]   (range: 1275348465.59)
%
% Effective contrast ratio (max range / min range): 85915.0618991
%
% USAGE:
%   f = transistor_9D(x)          % Evaluate function at point x (9D vector)
%   [lb, ub] = transistor_9D(n)   % Get bounds for dimension n (must be 9)
%   info = transistor_9D()        % Get complete problem information

nloc = 9;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [3;3;3;3;3;3;3;3;3];
lb_work = [-892907440.6673961;-325216597.0724246;-66650.01875769443;-11418.68794060073;-837085610.2400337;-382419663.414867;-59481.85819983898;-187265039.6239402;-981037281.2234801];
ub_work = [267872232.2002187;97564979.12172739;19995.00562730833;3425.606382180219;251125683.0720102;114725899.0244601;17844.55745995169;56179511.88718207;294311184.3670441];
scale_factors = [89290744.0667396;32521659.70724246;6665.001875769443;1141.868794060073;83708561.02400337;38241966.3414867;5948.185819983898;18726503.96239402;98103728.12234801];
contrast_ratio = 85915.0618991229;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'transistor';
    info.source_p = 62;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
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
    if arg1 ~= nloc, error('This instance is 9D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 9 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = transistor_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = transistor_orig(x);
%
% Purpose:
%
%    Function transistor is the function described in
%    Montaz et al. (2005).
%
%    dim = 9
%    Minimum global value = 0
%    Global minima near [0.9 0.45 1 2 8 8 5 1 2]'
%    Local minima
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
g = [0.485 0.752 0.869 0.982
    0.369 1.254 0.703 1.455
    5.2095 10.0677 22.9274 20.2153
    23.3037 101.779 111.461 191.267
    28.5132 111.8467 134.3884 211.4823]';
alfa = (1-x(1)*x(2))*x(3).*(exp(x(5).*(g(:,1)-g(:,3).*x(7)*0.001-...
       g(:,5).*x(8)*0.001))-1)-g(:,5)+g(:,4).*x(2);
beta = (1-x(1)*x(2))*x(4).*(exp(x(6).*(g(:,1)-g(:,2)-g(:,3).*x(7)*0.001+...
       g(:,4).*x(9)*0.001))-1)-g(:,5).*x(1)+g(:,4);
f = (x(1)*x(3)-x(2)*x(4))^2 + sum(alfa.^2+beta.^2);
%
% End of transistor.

