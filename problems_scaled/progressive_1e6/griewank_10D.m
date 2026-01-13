function varargout = griewank_10D(varargin)
%GRIEWANK_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   griewank (source instance p=22)
% Dimension: n = 10
% Strategy folder: progressive (kappa = 1000000)
% Original bound tag: bound(p) = 400
% Effective contrast: 1000000
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = griewank_orig(x_orig)

nloc = 10;
lb_orig = [-400;-400;-400;-400;-400;-400;-400;-400;-400;-400];
ub_orig = [400;400;400;400;400;400;400;400;400;400];
lb_work = [-400;-4000;-40000;-400000;-40;-4;-0.4;-400;-4000;-40000];
ub_work = [400;4000;40000;400000;40;4;0.4;400;4000;40000];
scale_factors = [1;10;100;1000;0.1;0.01;0.001;1;10;100];
contrast_ratio = 1000000;

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 22;
    info.dimension = nloc;
    info.strategy = 'progressive';
    info.kappa = 1000000;
    info.bound_p = 400;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Griewank (n=10): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
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
varargout{1} = griewank_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = griewank_orig(x);
%
% Purpose:
%
%    Function griewank is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -600 <= x <= 600 (Huyer and Neumaier (1999))
%                   -400 <= x <= 400 (Storn and Price (1997))
%    Cases considered: n = 10 Storn and Price (1997)
%                      n = 5 Huyer and Neumaier (1999), ICEO
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
n = size(x,1);
f = 1 + 1/4000*sum(x.^2,1)- prod(cos(x./sqrt([1:n]')),1);
%
% End of griewank.

