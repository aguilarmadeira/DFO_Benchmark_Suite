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
% Strategy folder: sobol_digit_oscillatory (kappa = 100000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 96771.82632112667
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = paviani_orig(x_orig)

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001];
lb_work = [19528.73664784381;62039343.10982545;12180.30158566732;3547.148690844169;15458.70140520449;90124789.01007479;131133252.3348034;1355.076754464178;18642.32793984628;68839610.97113925];
ub_work = [97585.12630774126;310010690.5323063;60864.9852849013;17725.10732621232;77247.15409827074;450353705.8029676;655273058.5185905;6771.32057365683;93155.74066492902;343991639.2305955];
scale_factors = [9759.488579632089;31004169.47017765;6087.10723921405;1772.688001421374;7725.487958622935;45039874.56775353;65533859.23778283;677.1997773434173;9316.505717064609;34402604.1834779];
contrast_ratio = 96771.82632112667;

if nargin == 0
    info.name = mfilename;
    info.problem = 'paviani';
    info.source_p = 36;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -45.778;
    info.x_global_min_orig = [9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001];
    info.x_global_min_work = [91260.97770813968;289919988.7156312;56920.5397938906;16576.40550129127;72241.03790108308;421167867.0830634;612807117.7325072;6332.495117938296;87118.64496027117;321698751.7197019];
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

