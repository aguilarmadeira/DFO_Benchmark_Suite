function varargout = paviani_10D(varargin)
%PAVIANI_10D  paviani 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [298727.414062, 1074836.46328]   (range: 776109.049219)
%   x2   ∈ [1298.22691406, 4671.08661328]   (range: 3372.85969922)
%   x3   ∈ [798477.164062, 2872961.53828]   (range: 2074484.37422)
%   x4   ∈ [1797.97666406, 6469.21168828]   (range: 4671.23502422)
%   x5   ∈ [173789.976562, 625305.194531]   (range: 451515.217969)
%   x6   ∈ [1173.28947656, 4221.55534453]   (range: 3048.26586797)
%   x7   ∈ [673539.726562, 2423430.26953]   (range: 1749890.54297)
%   x8   ∈ [1673.03922656, 6019.68041953]   (range: 4346.64119297)
%   x9   ∈ [423664.851562, 1524367.73203]   (range: 1100702.88047)
%   x10  ∈ [1423.16435156, 5120.61788203]   (range: 3697.45353047)
%
% Effective contrast ratio (max range / min range): 680.545747672
%
% USAGE:
%   f = paviani_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = paviani_10D(n)   % Get bounds for dimension n (must be 10)
%   info = paviani_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997];
lb_work = [298727.4140624999;1298.226914062499;798477.1640624998;1797.9766640625;173789.9765624999;1173.2894765625;673539.7265624999;1673.0392265625;423664.8515624999;1423.1643515625];
ub_work = [1074836.46328125;4671.08661328125;2872961.53828125;6469.21168828125;625305.19453125;4221.55534453125;2423430.26953125;6019.68041953125;1524367.73203125;5120.61788203125];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

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
    info.global_min_known = false;
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

