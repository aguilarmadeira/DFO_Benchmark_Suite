function varargout = langerman_10D(varargin)
%LANGERMAN_10D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   langerman (source instance p=29)
% Dimension: n = 10
% Strategy folder: extreme (kappa = 100000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 100000000
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = langerman_orig(x_orig)

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [10;10;10;10;10;10;10;10;10;10];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [10;10;10;10;10;1000000000;1000000000;1000000000;1000000000;1000000000];
scale_factors = [1;1;1;1;1;100000000;100000000;100000000;100000000;100000000];
contrast_ratio = 100000000;

if nargin == 0
    info.name = mfilename;
    info.problem = 'langerman';
    info.source_p = 29;
    info.dimension = nloc;
    info.strategy = 'extreme';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -0.965;
    info.x_global_min_orig = [8.074;8.776999999999999;3.467;1.867;6.708;6.349;4.534;0.276;7.633;1.567];
    info.x_global_min_work = [8.074;8.776999999999999;3.467;1.867;6.708;634900000;453399999.9999999;27600000;763300000;156700000];
    info.global_min_note = 'Langerman (10D): f*=-0.965. Ref: Ali et al. (2005).';
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
varargout{1} = langerman_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = langerman_orig(x);
%
% Purpose:
%
%    Function langerman is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -0.965, for n = 5, 10
%    Global minima: [8.074 8.777 3.467 1.867 6.708]' for n = 5
%                   [8.074 8.777 3.467 1.867 6.708 6.349 4.534...
%                    0.276 7.633 1.567]' for n = 10
%    Local minima 
%    Search domain: 0 <= x <= 10 (Montaz et al. (2005))
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
dim = size(x,1);
A   = [0.806 9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020
0.517 9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374
0.100 8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982
0.908 2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426
0.965 8.074 8.777 3.467 1.867 6.708 6.349 4.534 0.276 7.633 1.567];
d = sum((repmat(x',5,1)-A(:,2:11)).^2,2);
f = -sum(A(:,1).*cos(pi.*d).*exp(-d./pi),1);
%
% End of langerman.

