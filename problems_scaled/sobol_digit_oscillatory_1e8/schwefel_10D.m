function varargout = schwefel_10D(varargin)
%SCHWEFEL_10D  schwefel 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-147739.956072, 147739.956072]   (range: 295479.912144)
%   x2   ∈ [-2594636.49425, 2594636.49425]   (range: 5189272.9885)
%   x3   ∈ [-12566149640.9, 12566149640.9]   (range: 25132299281.7)
%   x4   ∈ [-38098260376.2, 38098260376.2]   (range: 76196520752.3)
%   x5   ∈ [-6316848093.42, 6316848093.42]   (range: 12633696186.8)
%   x6   ∈ [-31848982664.6, 31848982664.6]   (range: 63697965329.3)
%   x7   ∈ [-2022561.73856, 2022561.73856]   (range: 4045123.47712)
%   x8   ∈ [-44694606597.4, 44694606597.4]   (range: 89389213194.9)
%   x9   ∈ [-4630709205.56, 4630709205.56]   (range: 9261418411.13)
%   x10  ∈ [-29092022083.5, 29092022083.5]   (range: 58184044167 )
%
% Effective contrast ratio (max range / min range): 302522.132711
%
% Known global minimum (WORK-space):
%   x* = [124387.7944913293;2184521.503912171;10579911356.63456;32076350285.63306;5318390659.972579;26814849657.3202;1702870.371502045;37630060872.67728;3898767268.688183;24493661433.73133]
%   f* = -4189.829
%
% USAGE:
%   f = schwefel_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = schwefel_10D(n)   % Get bounds for dimension n (must be 10)
%   info = schwefel_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-500;-500;-500;-500;-500;-500;-500;-500;-500;-500];
ub_orig = [500;500;500;500;500;500;500;500;500;500];
lb_work = [-147739.956071947;-2594636.494247875;-12566149640.85757;-38098260376.16699;-6316848093.42426;-31848982664.64965;-2022561.738559239;-44694606597.44688;-4630709205.563481;-29092022083.50803];
ub_work = [147739.956071947;2594636.494247875;12566149640.85757;38098260376.16699;6316848093.42426;31848982664.64965;2022561.738559239;44694606597.44688;4630709205.563481;29092022083.50803];
scale_factors = [295.479912143894;5189.27298849575;25132299.28171514;76196520.75233398;12633696.18684852;63697965.32929931;4045.123477118478;89389213.19489376;9261418.411126962;58184044.16701605];
contrast_ratio = 302522.1327105399;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'schwefel';
    info.source_p = 46;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 500;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -4189.829;
    info.x_global_min_orig = [420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687;420.9687];
    info.x_global_min_work = [124387.7944913293;2184521.503912171;10579911356.63456;32076350285.63306;5318390659.972579;26814849657.3202;1702870.371502045;37630060872.67728;3898767268.688183;24493661433.73133];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Schwefel (n=10): x*=420.9687, f*=-418.9829*n. Ref: Ali et al. (2005).';
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
varargout{1} = schwefel_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = schwefel_orig(x);
%
% Purpose:
%
%    Function schwefel is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -418.9829*n
%    Global minima = 420.97*ones(n,1)
%    Local minima 
%    Search domain: -500 <= x <= 500 (Montaz et al. (2005))
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
f = -sum(x.*sin(sqrt(abs(x))),1);
%
% End of schwefel.

