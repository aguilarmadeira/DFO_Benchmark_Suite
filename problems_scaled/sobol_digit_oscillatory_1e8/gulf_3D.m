function varargout = gulf_3D(varargin)
%GULF_3D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   gulf (source instance p=24)
% Dimension: n = 3
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 15941.66048462197
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = gulf_orig(x_orig)

nloc = 3;
lb_orig = [0.1;0;0];
ub_orig = [100;25.6;5];
lb_work = [6361335.658449173;0;0];
ub_work = [6361335658.449685;102153.8459016891;44016.61782978238];
scale_factors = [63613356.58449685;3990.38460553473;8803.323565956476];
contrast_ratio = 15941.66048462197;

if nargin == 0
    info.name = mfilename;
    info.problem = 'gulf';
    info.source_p = 24;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [50;25;1.5];
    info.x_global_min_work = [3180667829.224842;99759.61513836826;13204.98534893471];
    info.global_min_note = 'Gulf (3D): x*=(50,25,1.5), f*=0. Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 3D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 3 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = gulf_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = gulf_orig(x);
%
% Purpose:
%
%    Function gulf is the function described in
%    Montaz et al. (2005).
%
%    dim = 3
%    Minimum global value = 0
%    Global minima = [50 25 1.5]'
%    Local minima 
%    Search domain: 0.1 <= x(1) <= 100 
%                   0 <= x(2) <= 25.6
%                   0 <= x(3) <= 5    (Montaz et al. (2005))
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
u = 25+(-50*log(0.01.*[1:99])).^(2/3);
f = sum((exp(-(u-x(2)).^x(3)./x(1))-0.01.*[1:99]).^2,2);
%
% End of gulf.

