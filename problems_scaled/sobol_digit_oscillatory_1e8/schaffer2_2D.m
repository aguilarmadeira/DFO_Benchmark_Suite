function varargout = schaffer2_2D(varargin)
%SCHAFFER2_2D  schaffer2 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-714157.331641, 714157.331641]   (range: 1428314.66328)
%   x2   ∈ [-200932.689434, 200932.689434]   (range: 401865.378869)
%
% Effective contrast ratio (max range / min range): 3.55421177933
%
% Known global minimum (WORK-space):
%   x* = [0;0]
%   f* = 0
%
% USAGE:
%   f = schaffer2_2D(x)          % Evaluate function at point x (2D vector)
%   [lb, ub] = schaffer2_2D(n)   % Get bounds for dimension n (must be 2)
%   info = schaffer2_2D()        % Get complete problem information

nloc = 2;
lb_orig = [-100;-100];
ub_orig = [100;100];
lb_work = [-714157.3316411743;-200932.6894344212];
ub_work = [714157.3316411743;200932.6894344212];
scale_factors = [7141.573316411743;2009.326894344213];
contrast_ratio = 3.554211779334467;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'schaffer2';
    info.source_p = 45;
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
    info.x_global_min_orig = [0;0];
    info.x_global_min_work = [0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Schaffer2 (2D): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = schaffer2_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schaffer2_orig(x);
%
% Purpose:
%
%    Function schaffer2 is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
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
f = (x(1)^2+x(2)^2)^0.25*(sin(50*(x(1)^2+x(2)^2)^0.1)^2+1);
%
% End of schaffer2.

