function varargout = fifteenn_local_minima_10D(varargin)
%FIFTEENN_LOCAL_MINIMA_10D  fifteenn_local_minima 10D test problem (heterogeneous WORK-space wrapper).
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
%   x10  ∈ [-3131.875   , 3131.875    ]   (range: 6263.75     )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [500500;250.75;750250.0000000009;125.875;625375.0000000009;375.625;875125;63.4375;562937.5000000009;313.1875000000005]
%   f* = 0
%
% USAGE:
%   f = fifteenn_local_minima_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = fifteenn_local_minima_10D(n)   % Get bounds for dimension n (must be 10)
%   info = fifteenn_local_minima_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-10;-10;-10;-10;-10;-10;-10;-10;-10;-10];
ub_orig = [10;10;10;10;10;10;10;10;10;10];
lb_work = [-5005000;-2507.5;-7502500;-1258.75;-6253750;-3756.25;-8751250;-634.375;-5629375;-3131.875];
ub_work = [5005000;2507.5;7502500;1258.75;6253750;3756.25;8751250;634.375;5629375;3131.875];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5;313.1875];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'fifteenn_local_minima';
    info.source_p = 19;
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
    info.x_global_min_orig = [1;1;1;1;1;1;1;1;1;1];
    info.x_global_min_work = [500500;250.75;750250.0000000009;125.875;625375.0000000009;375.625;875125;63.4375;562937.5000000009;313.1875000000005];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Fifteen Local Minima (n=10): x*=1, f*=0. Ref: Brachetti et al. (1997).';
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
varargout{1} = fifteenn_local_minima_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = fifteenn_local_minima_orig(x);
%
% Purpose:
%
%    Function fifteenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (15^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
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
f   = 0.1*sum(((x(1:dim-1)-1).^2).*(1+10*sin(3*pi.*x(2:dim)).^2),1);
f   = f + 0.1*(sin(3*pi*x(1))^2+(x(dim)-1)^2*(1+sin(2*pi*x(dim))^2));
%
% End of fifteenn_local_minima.

