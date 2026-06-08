function varargout = paviani_10D(varargin)
%PAVIANI_10D  paviani 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [1001500.5   , 3603449.85  ]   (range: 2601949.35  )
%   x2   ∈ [501.75075   , 1805.324775 ]   (range: 1303.574025 )
%   x3   ∈ [1501250.25  , 5401574.925 ]   (range: 3900324.675 )
%   x4   ∈ [251.875875  , 906.2622375 ]   (range: 654.3863625 )
%   x5   ∈ [1251375.375 , 4502512.3875]   (range: 3251137.0125)
%   x6   ∈ [751.625625  , 2704.3873125]   (range: 1952.7616875)
%   x7   ∈ [1751125.125 , 6300637.4625]   (range: 4549512.3375)
%   x8   ∈ [126.9384375 , 456.73096875]   (range: 329.79253125)
%   x9   ∈ [1126437.9375, 4052981.11875]   (range: 2926543.18125)
%   x10  ∈ [626.6881875 , 2254.85604375]   (range: 1628.16785625)
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% USAGE:
%   f = paviani_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = paviani_10D(n)   % Get bounds for dimension n (must be 10)
%   info = paviani_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997;7.1997];
lb_work = [1001500.5;501.7507499999998;1501250.25;251.875875;1251375.375;751.6256249999998;1751125.125;126.9384375;1126437.9375;626.6881874999999];
ub_work = [3603449.85;1805.324775;5401574.925;906.2622375000001;4502512.387499999;2704.3873125;6300637.4625;456.73096875;4052981.11875;2254.85604375];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5;313.1875];
contrast_ratio = 13795.07389162562;

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
    info.strategy = 'halton_oscillatory';
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

