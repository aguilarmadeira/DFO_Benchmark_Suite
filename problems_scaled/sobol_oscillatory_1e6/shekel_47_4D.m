function varargout = shekel_47_4D(varargin)
%SHEKEL_47_4D  Self-contained scaled test function.
%
% Wrapper/scaling formulation:
%   J. F. A. Madeira (2026)
%
% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   Journal of Global Optimization, 2026.
%
% Problem:   shekel_47 (source instance p=48)
% Dimension: n = 4
% Strategy folder: sobol_oscillatory (kappa = 1000000)
% Original bound tag: bound(p) = 0
% Effective contrast: 615.0520801974833
%
% Domain (scaled variables): x in [lb_work, ub_work] (see constants below)
% Mapping (as in create_scaled_wrapper.m):
%   t      = clip01((x - lb_work)./(ub_work - lb_work))
%   x_orig = lb_orig + t.*(ub_orig - lb_orig)
%   f      = shekel_47_orig(x_orig)

nloc = 4;
lb_orig = [0;0;0;0];
ub_orig = [10;10;10;10];
lb_work = [0;0;0;0];
ub_work = [1492890.625;6487.890625;3990390.625;8985.390625];
scale_factors = [149289.0625;648.7890625;399039.0625;898.5390625];
contrast_ratio = 615.0520801974833;

if nargin == 0
    info.name = mfilename;
    info.problem = 'shekel_47';
    info.source_p = 48;
    info.dimension = nloc;
    info.strategy = 'sobol_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -10.4029;
    info.x_global_min_orig = [4;4;4;4];
    info.x_global_min_work = [597156.25;2595.15625;1596156.25;3594.15625];
    info.global_min_note = 'Shekel-4,7 (4D): x*=(4,4,4,4), f*=-10.4029. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 4D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 4 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = shekel_47_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = shekel_47_orig(x);
%
% Purpose:
%
%    Function shekel_4 is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value  
%    Global minima = 4*ones(4,1)
%    Local minima (m local minima at x = A(:,1:4))
%    Search domain: 0 <= x <= 10 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
%    Cases considered: m = 5 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 7 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 10 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
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
m = 7;
A = [4 4 4 4 0.1
    1 1 1 1 0.2
    8 8 8 8 0.2
    6 6 6 6 0.4
    3 7 3 7 0.4
    2 9 2 9 0.6
    5 5 3 3 0.3
    8 1 8 1 0.7
    6 2 6 2 0.5
    7 3.6 7 3.6 0.5];
f = -sum(1./(sum((repmat(x',m,1)-A(1:m,1:4)).^2,2)+A(1:m,5)),1);
%
% End of shekel_4.

