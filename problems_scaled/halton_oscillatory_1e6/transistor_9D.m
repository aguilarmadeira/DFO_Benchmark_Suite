function varargout = transistor_9D(varargin)
%TRANSISTOR_9D  transistor 9D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-5005000    , 5005000     ]   (range: 10010000    )
%   x2   ∈ [-2507.5     , 2507.5      ]   (range: 5015        )
%   x3   ∈ [-7502500    , 7502500     ]   (range: 15005000    )
%   x4   ∈ [-1258.75    , 1258.75     ]   (range: 2517.5      )
%   x5   ∈ [-6253750    , 6253750     ]   (range: 12507500    )
%   x6   ∈ [-3756.25    , 3756.25     ]   (range: 7512.5      )
%   x7   ∈ [-8751250    , 8751250     ]   (range: 17502500    )
%   x8   ∈ [-634.375    , 634.375     ]   (range: 1268.75     )
%   x9   ∈ [-5629375    , 5629375     ]   (range: 11258750    )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [450450;112.8374999999996;750250.0000000009;251.75;5003000;3005;4375625;63.4375;1125875]
%   f* = 0
%
% USAGE:
%   f = transistor_9D(x)          % Evaluate function at point x (9D vector)
%   [lb, ub] = transistor_9D(n)   % Get bounds for dimension n (must be 9)
%   info = transistor_9D()        % Get complete problem information

nloc = 9;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10;10];
lb_work = [-5005000;-2507.5;-7502500;-1258.75;-6253750;-3756.25;-8751250;-634.375;-5629375];
ub_work = [5005000;2507.5;7502500;1258.75;6253750;3756.25;8751250;634.375;5629375];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'transistor';
    info.source_p = 62;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0.9;0.45;1;2;8;8;5;1;2];
    info.x_global_min_work = [450450;112.8374999999996;750250.0000000009;251.75;5003000;3005;4375625;63.4375;1125875];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Transistor (9D): x* approx, f*=0 (target). Ref: Ali et al. (2005).';
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

