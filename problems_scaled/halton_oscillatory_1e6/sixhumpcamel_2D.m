function varargout = sixhumpcamel_2D(varargin)
%SIXHUMPCAMEL_2D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   sixhumpcamel (source instance p=54)
% Dimension: n = 2
% Strategy folder: halton_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 1996.011964107677
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = sixhumpcamel_orig(x_orig)

nloc = 2;
lb_orig = [-3;-2];
ub_orig = [3;2];
lb_work = [-1501500;-501.5];
ub_work = [1501500;501.5];
scale_factors = [500500;250.75];
contrast_ratio = 1996.011964107677;

if nargin == 0
    info.name = mfilename;
    info.problem = 'sixhumpcamel';
    info.source_p = 54;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -1.031628453489877;
    info.x_global_min_orig = [0.08984200000000001 -0.08984200000000001;-0.712656 0.712656];
    info.x_global_min_work = [44965.92100000009 -44965.92100000009;-178.698492 178.698492];
    info.global_min_note = 'Six-Hump Camel (2D): 2 global minima, f*=-1.0316. Ref: Brachetti et al. (1997).';
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
varargout{1} = sixhumpcamel_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = sixhumpcamel_orig(x);
%
% Purpose:
%
%    Function sixhumpcamel is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = -1.0316285
%    Global minima: [-0.089842 0.712656]', [0.089842 -0.712656]'
%    Local minima (six)
%
%    Search domain: -2.5 <= x(1) <= 2.5 
%                   -1.5 <= x(2) <= 1.5 (Brachetti et al. (1997))
%
%    or
%
%                   -3 <= x(1) <= 3 
%                   -2 <= x(2) <= 2 (Huyer and Neumaier (1999))
%                                   (Huyer and Neumaier (2008))
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
f = (4 - 2.1*x(1)^2 + x(1)^4/3)*x(1)^2 + x(1)*x(2) + (4*x(2)^2 - 4)*x(2)^2;
%
% End of sixhumpcamel.

