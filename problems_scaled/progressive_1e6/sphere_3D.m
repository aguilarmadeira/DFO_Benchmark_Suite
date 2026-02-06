function varargout = sphere_3D(varargin)
%SPHERE_3D  sphere 3D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (PROGRESSIVE HETEROGENEITY):
%
%   x1   ∈ [-5.12       , 5.12        ]   (range: 10.24       )
%   x2   ∈ [-51.2       , 51.2        ]   (range: 102.4       )
%   x3   ∈ [-512        , 512         ]   (range: 1024        )
%
% Effective contrast ratio (max range / min range): 100
%
% Known global minimum (WORK-space):
%   x* = [0;0;0]
%   f* = 0
%
% USAGE:
%   f = sphere_3D(x)          % Evaluate function at point x (3D vector)
%   [lb, ub] = sphere_3D(n)   % Get bounds for dimension n (must be 3)
%   info = sphere_3D()        % Get complete problem information

nloc = 3;
lb_orig = [-5.12;-5.12;-5.12];
ub_orig = [5.12;5.12;5.12];
lb_work = [-5.12;-51.2;-512];
ub_work = [5.12;51.2;512];
scale_factors = [1;10;100];
contrast_ratio = 100;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'sphere';
    info.source_p = 55;
    info.dimension = nloc;
    info.strategy = 'progressive';
    info.kappa = 1000000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0];
    info.x_global_min_work = [0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Sphere (n=3): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
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
varargout{1} = sphere_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = sphere_orig(x);
%
% Purpose:
%
%    Function sphere is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (equal to global minima)
%    Search domain: -5.12 <= x <= 5.12 (Huyer and Neumaier (1999))
%    Cases considered: n = 3 Huyer and Neumaier (1999)
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
f = sum(x.^2,1);
%
% End of sphere.

