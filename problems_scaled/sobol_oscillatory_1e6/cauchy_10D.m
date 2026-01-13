function varargout = cauchy_10D(varargin)
%CAUCHY_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   cauchy (source instance p=7)
% Dimension: n = 10
% Strategy folder: sobol_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 680.5457476716454
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = cauchy_orig(x_orig)

nloc = 10;
lb_orig = [2;2;2;2;2;2;2;2;2;2];
ub_orig = [26;26;26;26;26;26;26;26;26;26];
lb_work = [298578.125;1297.578125;798078.125;1797.078125;173703.125;1172.703125;673203.125;1672.203125;423453.125;1422.453125];
ub_work = [3881515.625;16868.515625;10375015.625;23362.015625;2258140.625;15245.140625;8751640.625;21738.640625;5504890.625;18491.890625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625;86851.5625;586.3515625;336601.5625;836.1015625;211726.5625;711.2265625];
contrast_ratio = 680.5457476716454;

if nargin == 0
    info.name = mfilename;
    info.problem = 'cauchy';
    info.source_p = 7;
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
varargout{1} = cauchy_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = cauchy_orig(x);
%
% Purpose:
%
%    Function cauchy is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value
%    Global minima
%    Local minima
%    Search domain: y(1) <= x <= y(dim) (Brachetti et al. (1997))
%    Cases considered: n = 4 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
%                      n = 25 Brachetti et al. (1997)
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
if n == 4
    y = [3 7 12 17]';
else
    if n == 10
        y = [2 5 7 8 11 15 17 21 23 26]';
    else
        if n == 25
            y = [4.1 7.7 17.5 31.4 32.7 92.4 115.3 118.3 119 129.6 198.6];
            y = [y, 200.7 242.5 255 274.7 274.7 303.8 334.1 430 489.1 703.4 ];
            y = [y, 978 1656 1697.8 2745.6]';
        end
    end
end
f = -sum(log(pi)+log(1+(y-x).^2),1);
%
% End of cauchy.

