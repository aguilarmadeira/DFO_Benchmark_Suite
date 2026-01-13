function varargout = transistor_9D(varargin)
%TRANSISTOR_9D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   transistor (source instance p=62)
% Dimension: n = 9
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 10
% Effective contrast: 87073.14989944651
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = transistor_orig(x_orig)

nloc = 9;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10;10];
lb_work = [-236271589.5604044;-611394312.8170722;-46367.6628434852;-837319396.9476523;-11039.32285750715;-73748.96095573898;-33948.75326030552;-961228613.9601064;-127476060.8738314];
ub_work = [236271589.5604044;611394312.8170722;46367.6628434852;837319396.9476523;11039.32285750715;73748.96095573898;33948.75326030552;961228613.9601064;127476060.8738314];
scale_factors = [23627158.95604044;61139431.28170722;4636.76628434852;83731939.69476524;1103.932285750715;7374.896095573898;3394.875326030552;96122861.39601064;12747606.08738314];
contrast_ratio = 87073.14989944651;

if nargin == 0
    info.name = mfilename;
    info.problem = 'transistor';
    info.source_p = 62;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0.9;0.45;1;2;8;8;5;1;2];
    info.x_global_min_work = [21264443.06043643;27512744.07676816;4636.766284348523;167463879.3895304;8831.458286005723;58999.16876459117;16974.37663015276;96122861.39601076;25495212.17476627];
    info.global_min_note = 'Transistor (9D): x* approx, f*=0 (target). Ref: Ali et al. (2005).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 9D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 9 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = transistor_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = transistor_orig(x);
%
% Purpose:
%
%    Function transistor is the function described in
%    Montaz et al. (2005).
%
%    dim = 9
%    Minimum global value = 0
%    Global minima near [0.9 0.45 1 2 8 8 5 1 2]'
%    Local minima
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
g = [0.485 0.752 0.869 0.982
    0.369 1.254 0.703 1.455
    5.2095 10.0677 22.9274 20.2153
    23.3037 101.779 111.461 191.267
    28.5132 111.8467 134.3884 211.4823]';
alfa = (1-x(1)*x(2))*x(3).*(exp(x(5).*(g(:,1)-g(:,3).*x(7)*0.001-...
       g(:,5).*x(8)*0.001))-1)-g(:,5)+g(:,4).*x(2);
beta = (1-x(1)*x(2))*x(4).*(exp(x(6).*(g(:,1)-g(:,2)-g(:,3).*x(7)*0.001+...
       g(:,4).*x(9)*0.001))-1)-g(:,5).*x(1)+g(:,4);
f = (x(1)*x(3)-x(2)*x(4))^2 + sum(alfa.^2+beta.^2);
%
% End of transistor.

