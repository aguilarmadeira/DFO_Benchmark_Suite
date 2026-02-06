function varargout = fifteenn_local_minima_8D(varargin)
%FIFTEENN_LOCAL_MINIMA_8D  fifteenn_local_minima 8D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-361510850.542, 361510850.542]   (range: 723021701.084)
%   x2   ∈ [-64356.3165242, 64356.3165242]   (range: 128712.633048)
%   x3   ∈ [-8621.38032292, 8621.38032292]   (range: 17242.7606458)
%   x4   ∈ [-93051.7244978, 93051.7244978]   (range: 186103.448996)
%   x5   ∈ [-380517108.59, 380517108.59]   (range: 761034217.179)
%   x6   ∈ [-59910.6837054, 59910.6837054]   (range: 119821.367411)
%   x7   ∈ [-17168.8329209, 17168.8329209]   (range: 34337.6658418)
%   x8   ∈ [-82752.6427459, 82752.6427459]   (range: 165505.285492)
%
% Effective contrast ratio (max range / min range): 44136.4484963
%
% Known global minimum (WORK-space):
%   x* = [36151085.0541929;6435.631652419186;862.1380322920122;9305.17244977846;38051710.85895663;5991.068370537585;1716.883292090668;8275.264274587593]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_8D(x)          % Evaluate function at point x (8D vector)
%   [lb, ub] = fifteenn_local_minima_8D(n)   % Get bounds for dimension n (must be 8)
%   info = fifteenn_local_minima_8D()        % Get complete problem information

nloc = 8;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10];
lb_work = [-361510850.5419288;-64356.31652419175;-8621.380322920109;-93051.72449778448;-380517108.5895663;-59910.68370537585;-17168.83292090667;-82752.64274587587];
ub_work = [361510850.5419288;64356.31652419175;8621.380322920109;93051.72449778448;380517108.5895663;59910.68370537585;17168.83292090667;82752.64274587587];
scale_factors = [36151085.05419288;6435.631652419175;862.138032292011;9305.172449778449;38051710.85895663;5991.068370537585;1716.883292090667;8275.264274587587];
contrast_ratio = 44136.44849629868;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 18;
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
    info.x_global_min_orig = [1;1;1;1;1;1;1;1];
    info.x_global_min_work = [36151085.0541929;6435.631652419186;862.1380322920122;9305.17244977846;38051710.85895663;5991.068370537585;1716.883292090668;8275.264274587593];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=8): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 8D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 8 components.');
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

