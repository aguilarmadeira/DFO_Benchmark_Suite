function varargout = fifteenn_local_minima_6D(varargin)
%FIFTEENN_LOCAL_MINIMA_6D  fifteenn_local_minima 6D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-44198.5771005, 44198.5771005]   (range: 88397.154201)
%   x2   ∈ [-68907.4187552, 68907.4187552]   (range: 137814.83751)
%   x3   ∈ [-13131.6356687, 13131.6356687]   (range: 26263.2713374)
%   x4   ∈ [-878417964.366, 878417964.366]   (range: 1756835928.73)
%   x5   ∈ [-309338221.913, 309338221.913]   (range: 618676443.825)
%   x6   ∈ [-55636.7534723, 55636.7534723]   (range: 111273.506945)
%
% Effective contrast ratio (max range / min range): 66893.263454
%
% Known global minimum (WORK-space):
%   x* = [4419.857710048564;6890.741875516818;1313.163566867552;87841796.4366014;30933822.19127464;5563.675347232725]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_6D(x)          % Evaluate function at point x (6D vector)
%   [lb, ub] = fifteenn_local_minima_6D(n)   % Get bounds for dimension n (must be 6)
%   info = fifteenn_local_minima_6D()        % Get complete problem information

nloc = 6;
lb_orig = [-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10];
lb_work = [-44198.57710048557;-68907.41875516815;-13131.63566867551;-878417964.3660138;-309338221.912746;-55636.75347232721];
ub_work = [44198.57710048557;68907.41875516815;13131.63566867551;878417964.3660138;309338221.912746;55636.75347232721];
scale_factors = [4419.857710048557;6890.741875516815;1313.163566867551;87841796.43660137;30933822.1912746;5563.675347232721];
contrast_ratio = 66893.26345395124;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 17;
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
    info.x_global_min_orig = [1;1;1;1;1;1];
    info.x_global_min_work = [4419.857710048564;6890.741875516818;1313.163566867552;87841796.4366014;30933822.19127464;5563.675347232725];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=6): x*=1, f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 6D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 6 components.');
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

