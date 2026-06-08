function varargout = fletcher_powell_3D(varargin)
%FLETCHER_POWELL_3D  fletcher_powell 3D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-95703.9270694, 28711.1781208]   (range: 124415.10519)
%   x2   ∈ [-42141.4694634, 12642.440839]   (range: 54783.9103024)
%   x3   ∈ [-53601.1513596, 16080.3454079]   (range: 69681.4967675)
%
% Effective contrast ratio (max range / min range): 2.27101542229
%
% Known global minimum (WORK-space):
%   x* = [9570.392706935789;0;0]
%   f* = 0
%
% USAGE:
%   f = fletcher_powell_3D(x)          % Evaluate function at point x (3D vector)
%   [lb, ub] = fletcher_powell_3D(n)   % Get bounds for dimension n (must be 3)
%   info = fletcher_powell_3D()        % Get complete problem information

nloc = 3;
lb_orig = [-10;-10;-10];
ub_orig = [3;3;3];
lb_work = [-95703.92706935786;-42141.46946339063;-53601.1513595821];
ub_work = [28711.17812080736;12642.44083901719;16080.34540787463];
scale_factors = [9570.392706935785;4214.146946339063;5360.11513595821];
contrast_ratio = 2.271015422290822;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fletcher_powell';
    info.source_p = 20;
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
    info.x_global_min_orig = [1;0;0];
    info.x_global_min_work = [9570.392706935789;0;0];
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

