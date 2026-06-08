function varargout = griewank_5D(varargin)
%GRIEWANK_5D  griewank 5D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-300300000  , 90090000    ]   (range: 390390000   )
%   x2   ∈ [-150450     , 45135       ]   (range: 195585      )
%   x3   ∈ [-450150000  , 135045000   ]   (range: 585195000   )
%   x4   ∈ [-75525      , 22657.5     ]   (range: 98182.5     )
%   x5   ∈ [-375225000  , 112567500   ]   (range: 487792500   )
%
% Effective contrast ratio (max range / min range): 5960.27805362
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = griewank_5D(x)          % Evaluate function at point x (5D vector)
%   [lb, ub] = griewank_5D(n)   % Get bounds for dimension n (must be 5)
%   info = griewank_5D()        % Get complete problem information

nloc = 5;
lb_orig = [-600;-600;-600;-600;-600];
ub_orig = [180;180;180;180;180];
lb_work = [-300300000;-150450;-450150000;-75525;-375225000];
ub_work = [90090000;45135;135045000;22657.5;112567500];
scale_factors = [500500;250.75;750250;125.875;625375];
contrast_ratio = 5960.278053624627;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 23;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 600;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Griewank (n=5): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 5D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 5 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = griewank_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = griewank_orig(x);
%
% Purpose:
%
%    Function griewank is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -600 <= x <= 600 (Huyer and Neumaier (1999))
%                   -400 <= x <= 400 (Storn and Price (1997))
%    Cases considered: n = 10 Storn and Price (1997)
%                      n = 5 Huyer and Neumaier (1999), ICEO
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
f = 1 + 1/4000*sum(x.^2,1)- prod(cos(x./sqrt([1:n]')),1);
%
% End of griewank.

