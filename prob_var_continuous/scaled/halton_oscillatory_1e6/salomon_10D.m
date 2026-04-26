function varargout = salomon_10D(varargin)
%SALOMON_10D  salomon 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-50050000   , 15015000    ]   (range: 65065000    )
%   x2   ∈ [-25075      , 7522.5      ]   (range: 32597.5     )
%   x3   ∈ [-75025000   , 22507500    ]   (range: 97532500    )
%   x4   ∈ [-12587.5    , 3776.25     ]   (range: 16363.75    )
%   x5   ∈ [-62537500   , 18761250    ]   (range: 81298750    )
%   x6   ∈ [-37562.5    , 11268.75    ]   (range: 48831.25    )
%   x7   ∈ [-87512500   , 26253750    ]   (range: 113766250   )
%   x8   ∈ [-6343.75    , 1903.125    ]   (range: 8246.875    )
%   x9   ∈ [-56293750   , 16888125    ]   (range: 73181875    )
%   x10  ∈ [-31318.75   , 9395.625    ]   (range: 40714.375   )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [0;0;0;0;0;0;0;0;0;0]
%   f* = 0
%
% USAGE:
%   f = salomon_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = salomon_10D(n)   % Get bounds for dimension n (must be 10)
%   info = salomon_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-100;-100;-100;-100;-100;-100;-100;-100;-100;-100];
ub_orig = [30;30;30;30;30;30;30;30;30;30];
lb_work = [-50050000;-25075;-75025000;-12587.5;-62537500;-37562.5;-87512500;-6343.75;-56293750;-31318.75];
ub_work = [15015000;7522.5;22507500;3776.25;18761250;11268.75;26253750;1903.125;16888125;9395.625];
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
    info.problem = 'salomon';
    info.source_p = 43;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 100;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;0;0;0;0;0;0;0;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Salomon (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = salomon_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x);
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
%
% End of salomon.

