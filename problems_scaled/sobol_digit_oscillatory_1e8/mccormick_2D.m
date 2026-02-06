function varargout = mccormick_2D(varargin)
%MCCORMICK_2D  mccormick 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-7236.86626082, 19298.3100288]   (range: 26535.1762897)
%   x2   ∈ [-201809555.014, 201809555.014]   (range: 403619110.027)
%
% Effective contrast ratio (max range / min range): 13943.1590788
%
% Known global minimum (WORK-space):
%   x* = [-2639.043896445044;-104066460.5353765]
%   f* = -1.9133
%
% USAGE:
%   f = mccormick_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = mccormick_2D(n)   % Get bounds for dimension n (must be 2)
%   info = mccormick_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-1.5;-3];
ub_orig = [4;3];
lb_work = [-7236.866260818219;-201809555.0136583];
ub_work = [19298.31002884858;201809555.0136583];
scale_factors = [4824.577507212146;67269851.67121945];
contrast_ratio = 13943.15907883319;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'mccormick';
    info.source_p = 30;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1.9133;
    info.x_global_min_orig = [-0.547;-1.547];
    info.x_global_min_work = [-2639.043896445044;-104066460.5353765];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: McCormick (2D): f*=-1.9133. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 2D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 2 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = mccormick_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = mccormick_orig(x);
%
% Purpose:
%
%    Function mccormick is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -1.9133
%    Global minima = [-0.547 -1.547]'
%    Local minima (two)
%    Search domain: -1.5 <= x(1) <= 4
%                     -3 <= x(2) <= 3 (Montaz et al. (2005))
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
f = sin(x(1)+x(2))+(x(1)-x(2))^2-3/2*x(1)+5/2*x(2)+1;
%
% End of mccormick.

