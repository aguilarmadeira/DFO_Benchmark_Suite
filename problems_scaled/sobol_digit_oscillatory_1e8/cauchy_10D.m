function varargout = cauchy_10D(varargin)
%CAUCHY_10D  cauchy 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [33137468.2623, 430787087.41]   (range: 397649619.148)
%   x2   ∈ [148625857.593, 1932136148.71]   (range: 1783510291.12)
%   x3   ∈ [6301.63133946, 81921.207413]   (range: 75619.5760736)
%   x4   ∈ [153903990.053, 2000751870.68]   (range: 1846847880.63)
%   x5   ∈ [21247971.5232, 276223629.802]   (range: 254975658.279)
%   x6   ∈ [112086512.002, 1457124656.03]   (range: 1345038144.03)
%   x7   ∈ [76517191.6841, 994723491.894]   (range: 918206300.21)
%   x8   ∈ [19195.6641593, 249543.634071]   (range: 230347.969912)
%   x9   ∈ [4913.4127349, 63874.3655537]   (range: 58960.9528188)
%   x10  ∈ [13322.9237664, 173198.008964]   (range: 159875.085197)
%
% Effective contrast ratio (max range / min range): 31323.2366903
%
% USAGE:
%   f = cauchy_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = cauchy_10D(n)   % Get bounds for dimension n (must be 10)
%   info = cauchy_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2;2;2;2;2;2;2;2;2;2];
ub_orig = [26;26;26;26;26;26;26;26;26;26];
lb_work = [33137468.26231721;148625857.593172;6301.631339462787;153903990.0526788;21247971.5232245;112086512.0023094;76517191.68412703;19195.66415931351;4913.412734902227;13322.92376643354];
ub_work = [430787087.4101239;1932136148.711236;81921.20741301627;2000751870.684827;276223629.8019185;1457124656.030023;994723491.8936512;249543.6340710757;63874.36555372897;173198.0089636361];
scale_factors = [16568734.13115861;74312928.79658601;3150.815669731394;76951995.02633949;10623985.76161225;56043256.00115472;38258595.84206351;9597.832079656757;2456.706367451115;6661.461883216771];
contrast_ratio = 31323.23669034115;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'cauchy';
    info.source_p = 7;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
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

