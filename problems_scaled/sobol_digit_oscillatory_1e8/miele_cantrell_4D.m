function varargout = miele_cantrell_4D(varargin)
%MIELE_CANTRELL_4D  miele_cantrell 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-78655734.713, 78655734.713]   (range: 157311469.426)
%   x2   ∈ [-88050.7351752, 88050.7351752]   (range: 176101.47035)
%   x3   ∈ [-272334162.514, 272334162.514]   (range: 544668325.028)
%   x4   ∈ [-690559339.634, 690559339.634]   (range: 1381118679.27)
%
% Effective contrast ratio (max range / min range): 7842.74359844
%
% Known global minimum (WORK-space):
%   x* = [0;8805.073517523182;27233416.25142169;69055933.96335518]
%   f* = 0
%
% USAGE:
%   f = miele_cantrell_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = miele_cantrell_4D(n)   % Get bounds for dimension n (must be 4)
%   info = miele_cantrell_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [10;10;10;10];
lb_work = [-78655734.71298145;-88050.73517523179;-272334162.5142167;-690559339.6335516];
ub_work = [78655734.71298145;88050.73517523179;272334162.5142167;690559339.6335516];
scale_factors = [7865573.471298145;8805.073517523178;27233416.25142167;69055933.96335515];
contrast_ratio = 7842.743598440758;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'miele_cantrell';
    info.source_p = 32;
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
    info.x_global_min_orig = [0;1;1;1];
    info.x_global_min_work = [0;8805.073517523182;27233416.25142169;69055933.96335518];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Miele-Cantrell (4D): x*=(0,1,1,1), f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = miele_cantrell_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = miele_cantrell_orig(x);
%
% Purpose:
%
%    Function miele_cantrell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [0 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = (exp(x(1))-x(2))^4 + 100*(x(2)-x(3))^6 + tan(x(3)-x(4))^4 + x(1)^8;
%
% End of miele_cantrell.

