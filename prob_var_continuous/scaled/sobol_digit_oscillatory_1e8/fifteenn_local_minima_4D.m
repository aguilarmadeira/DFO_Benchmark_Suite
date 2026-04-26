function varargout = fifteenn_local_minima_4D(varargin)
%FIFTEENN_LOCAL_MINIMA_4D  fifteenn_local_minima 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-961282543.101, 288384762.93]   (range: 1249667306.03)
%   x2   ∈ [-253110100.319, 75933030.0957]   (range: 329043130.415)
%   x3   ∈ [-62480.4054232, 18744.1216269]   (range: 81224.5270501)
%   x4   ∈ [-160728783.329, 48218634.9987]   (range: 208947418.328)
%
% Effective contrast ratio (max range / min range): 15385.3441986
%
% Known global minimum (WORK-space):
%   x* = [96128254.31012642;25311010.03190139;6248.040542316579;16072878.33291158]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = fifteenn_local_minima_4D(n)   % Get bounds for dimension n (must be 4)
%   info = fifteenn_local_minima_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [3;3;3;3];
lb_work = [-961282543.1012644;-253110100.3190137;-62480.40542316585;-160728783.3291159];
ub_work = [288384762.9303793;75933030.09570411;18744.12162694976;48218634.99873478];
scale_factors = [96128254.31012644;25311010.03190137;6248.040542316585;16072878.33291159];
contrast_ratio = 15385.34419856453;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 16;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;1;1;1];
    info.x_global_min_work = [96128254.31012642;25311010.03190139;6248.040542316579;16072878.33291158];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=4): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = fifteenn_local_minima_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = fifteenn_local_minima_orig(x);
%
% Purpose:
%
%    Function fifteenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (15^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
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
f   = 0.1*sum(((x(1:dim-1)-1).^2).*(1+10*sin(3*pi.*x(2:dim)).^2),1);
f   = f + 0.1*(sin(3*pi*x(1))^2+(x(dim)-1)^2*(1+sin(2*pi*x(dim))^2));
%
% End of fifteenn_local_minima.

