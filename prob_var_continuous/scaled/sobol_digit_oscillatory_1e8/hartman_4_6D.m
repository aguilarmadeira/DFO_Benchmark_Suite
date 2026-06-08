function varargout = hartman_4_6D(varargin)
%HARTMAN_4_6D  hartman_4 6D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 749760.368697]   (range: 749760.368697)
%   x2   ∈ [0           , 45713673.8009]   (range: 45713673.8009)
%   x3   ∈ [0           , 21198241.2996]   (range: 21198241.2996)
%   x4   ∈ [0           , 58301315.458]   (range: 58301315.458)
%   x5   ∈ [0           , 1032.613699 ]   (range: 1032.613699 )
%   x6   ∈ [0           , 40183147.575]   (range: 40183147.575)
%
% Effective contrast ratio (max range / min range): 56459.9477178
%
% USAGE:
%   f = hartman_4_6D(x)          % Evaluate function at point x (6D vector)
%   [lb, ub] = hartman_4_6D(n)   % Get bounds for dimension n (must be 6)
%   info = hartman_4_6D()        % Get complete problem information

nloc = 6;
lb_orig = [0;0;0;0;0;0];
ub_orig = [0.65;0.65;0.65;0.65;0.65;0.65];
lb_work = [0;0;0;0;0;0];
ub_work = [749760.3686966163;45713673.8009389;21198241.29959079;58301315.4579775;1032.61369899569;40183147.57502275];
scale_factors = [1153477.490302487;70328728.92452139;32612678.92244737;89694331.47381154;1588.636459993369;61820227.03849654];
contrast_ratio = 56459.94771779687;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'hartman_4';
    info.source_p = 26;
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
    if arg1 ~= nloc, error('This instance is 6D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 6 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = hartman_4_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = hartman_4_orig(x);
%
% Purpose:
%
%    Function hartman_4 is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value
%    Global minima (one)
%    Local minima (4 local minima at x = A(:,n+2:2n+1))
%    Search domain: 0 <= x <= 1 (Brachetti et al. (1997))
%                               (Huyer and Neumaier (1999))
%                               (Huyer and Neumaier (2008))
%    Cases considered: n = 3 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      n = 6 Brachetti et al. (1997)
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
n = size(x,1);
if n == 3
    A = [3 10 30 1 0.3689 0.1170 0.2673
        0.1 10 35 1.2 0.4699 0.4387 0.7470
        3 10 30 3 0.1091 0.8732 0.5547
        0.1 10 35 3.2 0.03815 0.5743 0.8828];
else
    if n== 6
        A = [10 3 17 3.5 1.7 8 1 0.1312 0.1696 0.5569 0.0124 0.8283 0.5886
            0.05 10 17 0.1 8 14 1.2 0.2329 0.4135 0.8307 0.3736 0.1004 0.9991
            3 3.5 1.7 10 17 8 3 0.2348 0.1451 0.3522 0.2883 0.3047 0.6650
            17 8 0.05 10 0.1 14 3.2 0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
    end
end
f = -sum(A(:,n+1).*exp(-sum(A(:,1:n).*(repmat(x',4,1)-A(:,n+2:2*n+1)).^2,2)),1);
%
% End of hartman_4.

