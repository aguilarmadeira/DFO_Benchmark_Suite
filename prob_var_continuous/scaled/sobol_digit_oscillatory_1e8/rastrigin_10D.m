function varargout = rastrigin_10D(varargin)
%RASTRIGIN_10D  rastrigin 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-281037639.381, 84311291.8144]   (range: 365348931.196)
%   x2   ∈ [-73823630.6275, 22147089.1882]   (range: 95970719.8157)
%   x3   ∈ [-47903.6121869, 14371.0836561]   (range: 62274.6958429)
%   x4   ∈ [-14383.5621627, 4315.06864882]   (range: 18698.6308116)
%   x5   ∈ [-347153834.348, 104146150.304]   (range: 451299984.652)
%   x6   ∈ [-1195.30900194, 358.592700583]   (range: 1553.90170253)
%   x7   ∈ [-41316.9902348, 12395.0970704]   (range: 53712.0873052)
%   x8   ∈ [-20595.613774, 6178.6841322]   (range: 26774.2979062)
%   x9   ∈ [-30387.7477301, 9116.32431902]   (range: 39504.0720491)
%   x10  ∈ [-12716.0912385, 3814.82737155]   (range: 16530.91861 )
%
% Effective contrast ratio (max range / min range): 290430.201549
%
% Known global minimum (WORK-space):
%   x* = [0;-1.490116119384766e-08;0;-1.818989403545856e-12;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = rastrigin_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = rastrigin_10D(n)   % Get bounds for dimension n (must be 10)
%   info = rastrigin_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12];
ub_orig = [1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536];
lb_work = [-281037639.3812032;-73823630.62747142;-47903.61218687722;-14383.56216273422;-347153834.3475208;-1195.309001944348;-41316.99023478182;-20595.61377398605;-30387.74773007598;-12716.09123849275];
ub_work = [84311291.81436099;22147089.18824143;14371.08365606317;4315.068648820267;104146150.3042563;358.5927005833046;12395.09707043455;6178.684132195816;9116.324319022795;3814.827371547826];
scale_factors = [54890163.94164125;14418677.85692801;9356.174255249458;2809.289484909028;67803483.27100016;233.4587894422556;8069.724655230824;4022.58081523165;5935.106978530464;2483.611570018116];
contrast_ratio = 290430.2015485731;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'rastrigin';
    info.source_p = 40;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;-1.490116119384766e-08;0;-1.818989403545856e-12;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Rastrigin (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = rastrigin_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = rastrigin_orig(x);
%
% Purpose:
%
%    Function rastrigin is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -5.12 <= x <= 5.12 (Montaz et al. (2005))
%    Cases considered: n = 10 Montaz et al. (2005)
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
dim = size(x,1);
f   = 10*dim + sum(x.^2-10.*cos(2*pi.*x),1);
%
% End of rastrigin.

