function varargout = neumaier2_4D(varargin)
%NEUMAIER2_4D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   neumaier2 (source instance p=34)
% Dimension: n = 4
% Strategy folder: spatial_thermal (kappa = 90000)
% Original bound tag: bound(p) = 0
% Effective contrast: 300
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = neumaier2_orig(x_orig)

nloc = 4;
lb_orig = [0;0;0;0];
ub_orig = [4;4;4;4];
lb_work = [0;0;0;0];
ub_work = [4;4;1200;1200];
scale_factors = [1;1;300;300];
contrast_ratio = 300;

if nargin == 0
    info.name = mfilename;
    info.problem = 'neumaier2';
    info.source_p = 34;
    info.dimension = nloc;
    info.strategy = 'spatial_thermal';
    info.kappa = 90000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [1;2;2;3];
    info.x_global_min_work = [1;2;600;900];
    info.global_min_note = 'Neumaier2 (4D): x*=(1,2,2,3), f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = neumaier2_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = neumaier2_orig(x);
%
% Purpose:
%
%    Function neumaier2 is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = [1 2 2 3]'
%    Local minima
%    Search domain: 0 <= x <= n   (Montaz et al. (2005))
%    Cases considered: n = 4, b = [8 18 44 114]' Montaz et al. (2005)
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
b   = [8 18 44 114]';
aux = repmat([1:4]',1,4);
f   = sum((b - sum(repmat(x',4,1).^aux,2)).^2,1);
%
% End of neumaier2.

