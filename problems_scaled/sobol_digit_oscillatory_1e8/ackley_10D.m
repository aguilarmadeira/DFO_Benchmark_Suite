function varargout = ackley_10D(varargin)
%ACKLEY_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   ackley (source instance p=1)
% Dimension: n = 10
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 30
% Effective contrast: 149515.9300171582
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = ackley_orig(x_orig)

nloc = 10;
lb_orig = [-30;-30;-30;-30;-30;-30;-30;-30;-30;-30];
ub_orig = [30;30;30;30;30;30;30;30;30;30];
lb_work = [-35974.38712086231;-240956.6263231844;-1485952257.170651;-2038854039.18265;-62843.63784577073;-2703201531.597895;-100792.5327610986;-1579911090.063973;-18079.68910929879;-258218.9494962456];
ub_work = [35974.38712086231;240956.6263231844;1485952257.170651;2038854039.18265;62843.63784577073;2703201531.597895;100792.5327610986;1579911090.063973;18079.68910929879;258218.9494962456];
scale_factors = [1199.146237362077;8031.887544106147;49531741.90568837;67961801.30608833;2094.787928192357;90106717.71992984;3359.751092036619;52663703.00213243;602.6563036432931;8607.298316541521];
contrast_ratio = 149515.9300171582;

if nargin == 0
    info.name = mfilename;
    info.problem = 'ackley';
    info.source_p = 1;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 30;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Ackley (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = ackley_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = ackley_orig(x);
%
% Purpose:
%
%    Function ackley is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -30 <= x <= 30 (Montaz et al. (2005))
%    Cases considered: n = 10 Montaz et al. (2005)
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
f = -20*exp(-0.02*sqrt(1/n*sum(x.^2,1)))-exp(1/n*sum(cos(2*pi.*x),1))+20+exp(1);
%
% End of ackley.

