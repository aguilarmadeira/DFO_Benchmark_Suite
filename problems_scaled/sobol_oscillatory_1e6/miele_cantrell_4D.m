function varargout = miele_cantrell_4D(varargin)
%MIELE_CANTRELL_4D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   miele_cantrell (source instance p=32)
% Dimension: n = 4
% Strategy folder: sobol_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 10
% Effective contrast: 615.0520801974833
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = miele_cantrell_orig(x_orig)

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [10;10;10;10];
lb_work = [-1492890.625;-6487.890625;-3990390.625;-8985.390625];
ub_work = [1492890.625;6487.890625;3990390.625;8985.390625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625];
contrast_ratio = 615.0520801974833;

if nargin == 0
    info.name = mfilename;
    info.problem = 'miele_cantrell';
    info.source_p = 32;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;1;1;1];
    info.x_global_min_work = [0;648.7890625000009;399039.0625;898.5390625];
    info.global_min_note = 'Miele-Cantrell (4D): x*=(0,1,1,1), f*=0. Ref: Brachetti et al. (1997).';
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

