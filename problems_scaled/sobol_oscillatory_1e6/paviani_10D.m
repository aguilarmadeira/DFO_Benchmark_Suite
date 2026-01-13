function varargout = paviani_10D(varargin)
%PAVIANI_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   paviani (source instance p=36)
% Dimension: n = 10
% Strategy folder: sobol_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 680.5457476716454
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = paviani_orig(x_orig)

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001];
lb_work = [298727.4140624999;1298.226914062499;798477.1640624998;1797.976664062499;173789.9765624999;1173.2894765625;673539.7265624998;1673.0392265625;423664.8515624999;1423.1643515625];
ub_work = [1492741.3359375;6487.241835937501;3989991.5859375;8984.492085937502;868428.7734375;5862.9292734375;3365679.0234375;8360.1795234375;2117053.8984375;7111.5543984375];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

if nargin == 0
    info.name = mfilename;
    info.problem = 'paviani';
    info.source_p = 36;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -45.778;
    info.x_global_min_orig = [9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001];
    info.x_global_min_work = [1396002.0234375;6066.826523437501;3731414.2734375;8402.238773437501;812148.9609375;5482.973460937501;3147561.2109375;7818.385710937501;1979855.0859375;6650.6795859375];
    info.global_min_note = 'Paviani (10D): x*=9.351, f*=-45.778. Ref: Ali et al. (2005).';
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
varargout{1} = paviani_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = paviani_orig(x);
%
% Purpose:
%
%    Function paviani is the function described in
%    Montaz et al. (2005).
%
%    dim = 10
%    Minimum global value = -45.778
%    Global minima = 9.351*ones(10,1)
%    Local minima
%    Search domain: 2.001 <= x <= 9.999
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
f = sum((log(x-2)).^2+(log(10-x)).^2,1)-prod(x,1)^0.2;
%
% End of paviani.

