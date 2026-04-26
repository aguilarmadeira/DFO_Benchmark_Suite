function varargout = ackley_10D(varargin)
%ACKLEY_10D  ackley 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-15015000   , 4504500     ]   (range: 19519500    )
%   x2   ∈ [-7522.5     , 2256.75     ]   (range: 9779.25     )
%   x3   ∈ [-22507500   , 6752250     ]   (range: 29259750    )
%   x4   ∈ [-3776.25    , 1132.875    ]   (range: 4909.125    )
%   x5   ∈ [-18761250   , 5628375     ]   (range: 24389625    )
%   x6   ∈ [-11268.75   , 3380.625    ]   (range: 14649.375   )
%   x7   ∈ [-26253750   , 7876125     ]   (range: 34129875    )
%   x8   ∈ [-1903.125   , 570.9375    ]   (range: 2474.0625   )
%   x9   ∈ [-16888125   , 5066437.5   ]   (range: 21954562.5  )
%   x10  ∈ [-9395.625   , 2818.6875   ]   (range: 12214.3125  )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = ackley_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = ackley_10D(n)   % Get bounds for dimension n (must be 10)
%   info = ackley_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-30;-30;-30;-30;-30;-30;-30;-30;-30;-30];
ub_orig = [9;9;9;9;9;9;9;9;9;9];
lb_work = [-15015000;-7522.5;-22507500;-3776.25;-18761250;-11268.75;-26253750;-1903.125;-16888125;-9395.625];
ub_work = [4504500;2256.75;6752250;1132.875;5628375;3380.625;7876125;570.9375;5066437.5;2818.6875];
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
    info.problem = 'ackley';
    info.source_p = 1;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 30;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Ackley (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = ackley_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = ackley_orig(x);
%
% Purpose:
%
%    Function ackley is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -30 <= x <= 30 (Montaz et al. (2005))
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
n = size(x,1);
f = -20*exp(-0.02*sqrt(1/n*sum(x.^2,1)))-exp(1/n*sum(cos(2*pi.*x),1))+20+exp(1);
%
% End of ackley.

