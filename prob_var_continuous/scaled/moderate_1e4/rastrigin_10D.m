function varargout = rastrigin_10D(varargin)
%RASTRIGIN_10D  rastrigin 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (MODERATE HETEROGENEITY):
%
%   x1   ∈ [-5.12       , 1.536       ]   (range: 6.656       )
%   x2   ∈ [-512        , 153.6       ]   (range: 665.6       )
%   x3   ∈ [-0.0512     , 0.01536     ]   (range: 0.06656     )
%   x4   ∈ [-5.12       , 1.536       ]   (range: 6.656       )
%   x5   ∈ [-512        , 153.6       ]   (range: 665.6       )
%   x6   ∈ [-0.0512     , 0.01536     ]   (range: 0.06656     )
%   x7   ∈ [-5.12       , 1.536       ]   (range: 6.656       )
%   x8   ∈ [-512        , 153.6       ]   (range: 665.6       )
%   x9   ∈ [-0.0512     , 0.01536     ]   (range: 0.06656     )
%   x10  ∈ [-5.12       , 1.536       ]   (range: 6.656       )
%
% Effective contrast ratio (max range / min range): 10000
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = rastrigin_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = rastrigin_10D(n)   % Get bounds for dimension n (must be 10)
%   info = rastrigin_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12];
ub_orig = [1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536];
lb_work = [-5.12;-512;-0.0512;-5.12;-512;-0.0512;-5.12;-512;-0.0512;-5.12];
ub_work = [1.536;153.6;0.01536000000000001;1.536;153.6;0.01536000000000001;1.536;153.6;0.01536000000000001;1.536];
scale_factors = [1;100;0.01;1;100;0.01;1;100;0.01;1];
contrast_ratio = 10000;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'rastrigin';
    info.source_p = 40;
    info.dimension = nloc;
    info.strategy = 'moderate';
    info.kappa = 10000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Rastrigin (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = rastrigin_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = rastrigin_orig(x);
%
% Purpose:
%
%    Function rastrigin is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -5.12 <= x <= 5.12 (Montaz et al. (2005))
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
f   = 10*dim + sum(x.^2-10.*cos(2*pi.*x),1);
%
% End of rastrigin.

