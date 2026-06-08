function varargout = sinusoidal_10D(varargin)
%SINUSOIDAL_10D  sinusoidal 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 620466.693819]   (range: 620466.693819)
%   x2   ∈ [0           , 143528867.263]   (range: 143528867.263)
%   x3   ∈ [0           , 924769.431014]   (range: 924769.431014)
%   x4   ∈ [0           , 3609334450.4]   (range: 3609334450.4)
%   x5   ∈ [0           , 8694518783.86]   (range: 8694518783.86)
%   x6   ∈ [0           , 267907.875421]   (range: 267907.875421)
%   x7   ∈ [0           , 1114346.46057]   (range: 1114346.46057)
%   x8   ∈ [0           , 545940.635616]   (range: 545940.635616)
%   x9   ∈ [0           , 699341.672455]   (range: 699341.672455)
%   x10  ∈ [0           , 135516.748696]   (range: 135516.748696)
%
% Effective contrast ratio (max range / min range): 64158.2599018
%
% USAGE:
%   f = sinusoidal_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = sinusoidal_10D(n)   % Get bounds for dimension n (must be 10)
%   info = sinusoidal_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [117;117;117;117;117;117;117;117;117;117];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [620466.6938188393;143528867.2632016;924769.431013884;3609334450.40402;8694518783.855522;267907.8754213403;1114346.46057429;545940.6356156491;699341.6724553917;135516.7486955722];
scale_factors = [5303.134135203754;1226742.45524104;7904.012230887898;30849012.39661555;74312126.35773951;2289.810901037097;9524.328722857177;4666.159278766231;5977.279251755484;1158.262809363865];
contrast_ratio = 64158.25990178585;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'sinusoidal';
    info.source_p = 53;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
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
varargout{1} = sinusoidal_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = sinusoidal_orig(x);
%
% Purpose:
%
%    Function sinusoidal is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -(A+1)
%    Global minima = (90+z)*ones(n,1)
%    Local minima (several)
%    Search domain: 0 <= x <= 180 (Montaz et al. (2005))
%    Cases considered: n = 10, 20
%                      A = 2.5;B = 5;z = 30;Montaz et al. (2005)
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
A = 2.5;
B = 5;
z = 30;
f = -(A*prod(sind(x-z),1)+prod(sind(B.*(x-z)),1));
%
% End of sinusoidal.

