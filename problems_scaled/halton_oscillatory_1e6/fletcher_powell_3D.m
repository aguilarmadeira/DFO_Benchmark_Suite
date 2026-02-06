function varargout = fletcher_powell_3D(varargin)
%FLETCHER_POWELL_3D  fletcher_powell 3D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-5005000    , 5005000     ]   (range: 10010000    )
%   x2   ∈ [-2507.5     , 2507.5      ]   (range: 5015        )
%   x3   ∈ [-7502500    , 7502500     ]   (range: 15005000    )
%
% Effective contrast ratio (max range / min range): 2992.02392822
%
% Known global minimum (WORK-space):
%   x* = [500500;0;0]
%   f* = 0
%
% USAGE:
%   f = fletcher_powell_3D(x)          % Evaluate function at point x (3D vector)
%   [lb, ub] = fletcher_powell_3D(n)   % Get bounds for dimension n (must be 3)
%   info = fletcher_powell_3D()        % Get complete problem information

nloc = 3;
lb_orig = [-10;-10;-10];
ub_orig = [10;10;10];
lb_work = [-5005000;-2507.5;-7502500];
ub_work = [5005000;2507.5;7502500];
scale_factors = [500500;250.75;750250];
contrast_ratio = 2992.023928215354;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fletcher_powell';
    info.source_p = 20;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;0;0];
    info.x_global_min_work = [500500;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fletcher-Powell (3D): x*=(1,0,0), f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 3D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 3 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = fletcher_powell_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = fletcher_powell_orig(x);
%
% Purpose:
%
%    Function fletcher_powell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 3
%    Minimum global value = 0
%    Global minima = [1 0 0]'
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
if x(1) == 0
    f = +inf;
else
   if x(1)>0
      teta = 1/(2*pi)*atan(x(2)/x(1));
   else
      teta = 1/(2*pi)*atan(x(2)/x(1))+0.5;
   end
   f = 100* ((x(3)-10*teta)^2+(norm(x,2)-1)^2)+x(3)^2;
end
%
% End of fletcher_powell.

