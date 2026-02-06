function varargout = fifteenn_local_minima_10D(varargin)
%FIFTEENN_LOCAL_MINIMA_10D  fifteenn_local_minima 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-46325.1963707, 46325.1963707]   (range: 92650.3927414)
%   x2   ∈ [-88407.236888, 88407.236888]   (range: 176814.473776)
%   x3   ∈ [-147341776.063, 147341776.063]   (range: 294683552.127)
%   x4   ∈ [-693260177.659, 693260177.659]   (range: 1386520355.32)
%   x5   ∈ [-258734378.633, 258734378.633]   (range: 517468757.265)
%   x6   ∈ [-839652479.733, 839652479.733]   (range: 1679304959.47)
%   x7   ∈ [-6676.2636572, 6676.2636572]   (range: 13352.5273144)
%   x8   ∈ [-522826858.392, 522826858.392]   (range: 1045653716.78)
%   x9   ∈ [-39049.6523544, 39049.6523544]   (range: 78099.3047088)
%   x10  ∈ [-965678005.315, 965678005.315]   (range: 1931356010.63)
%
% Effective contrast ratio (max range / min range): 144643.479482
%
% Known global minimum (WORK-space):
%   x* = [4632.519637068108;8840.723688800659;14734177.60633656;69326017.76590669;25873437.86326867;83965247.97325504;667.6263657196896;52282685.83915204;3904.965235441865;96567800.53147376]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = fifteenn_local_minima_10D(n)   % Get bounds for dimension n (must be 10)
%   info = fifteenn_local_minima_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10;10;10];
lb_work = [-46325.19637068104;-88407.23688800648;-147341776.0633654;-693260177.659066;-258734378.6326868;-839652479.7325494;-6676.263657196891;-522826858.3915203;-39049.65235441862;-965678005.3147368];
ub_work = [46325.19637068104;88407.23688800648;147341776.0633654;693260177.659066;258734378.6326868;839652479.7325494;6676.263657196891;522826858.3915203;39049.65235441862;965678005.3147368];
scale_factors = [4632.519637068104;8840.723688800648;14734177.60633654;69326017.7659066;25873437.86326868;83965247.97325495;667.6263657196891;52282685.83915203;3904.965235441862;96567800.53147368];
contrast_ratio = 144643.4794817837;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 19;
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
    info.x_global_min_orig = [1;1;1;1;1;1;1;1;1;1];
    info.x_global_min_work = [4632.519637068108;8840.723688800659;14734177.60633656;69326017.76590669;25873437.86326867;83965247.97325504;667.6263657196896;52282685.83915204;3904.965235441865;96567800.53147376];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=10): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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

