function varargout = langerman_10D(varargin)
%LANGERMAN_10D  langerman 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [0           , 585374337.681]   (range: 585374337.681)
%   x2   ∈ [0           , 40648.0859047]   (range: 40648.0859047)
%   x3   ∈ [0           , 787089452.878]   (range: 787089452.878)
%   x4   ∈ [0           , 20586.0514648]   (range: 20586.0514648)
%   x5   ∈ [0           , 695615164.732]   (range: 695615164.732)
%   x6   ∈ [0           , 31163.1052208]   (range: 31163.1052208)
%   x7   ∈ [0           , 93245.5835678]   (range: 93245.5835678)
%   x8   ∈ [0           , 75892775.5628]   (range: 75892775.5628)
%   x9   ∈ [0           , 523850454.041]   (range: 523850454.041)
%   x10  ∈ [0           , 48352.1383479]   (range: 48352.1383479)
%
% Effective contrast ratio (max range / min range): 38234.1146977
%
% Known global minimum (WORK-space):
%   x* = [472631240.243553;35676.82499854233;272883913.3128651;3843.415808478328;466618652.5020162;19785.45550466852;42277.54758962381;2094640.605532714;399855051.5693428;7576.780079111181]
%   f* = -0.965
%
% USAGE:
%   f = langerman_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = langerman_10D(n)   % Get bounds for dimension n (must be 10)
%   info = langerman_10D()        % Get complete problem information

nloc = 10;
lb_orig = [0;0;0;0;0;0;0;0;0;0];
ub_orig = [10;10;10;10;10;10;10;10;10;10];
lb_work = [0;0;0;0;0;0;0;0;0;0];
ub_work = [585374337.6808931;40648.08590468535;787089452.8781803;20586.0514648009;695615164.7316878;31163.10522077259;93245.58356776315;75892775.56277947;523850454.0408002;48352.1383478697];
scale_factors = [58537433.7680893;4064.808590468535;78708945.28781803;2058.60514648009;69561516.47316878;3116.310522077259;9324.558356776315;7589277.556277948;52385045.40408002;4835.21383478697];
contrast_ratio = 38234.11469771106;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'langerman';
    info.source_p = 29;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -0.965;
    info.x_global_min_orig = [8.074;8.776999999999999;3.467;1.867;6.708;6.349;4.534;0.276;7.633;1.567];
    info.x_global_min_work = [472631240.243553;35676.82499854233;272883913.3128651;3843.415808478328;466618652.5020162;19785.45550466852;42277.54758962381;2094640.605532714;399855051.5693428;7576.780079111181];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Langerman (10D): f*=-0.965. Ref: Ali et al. (2005).';
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
varargout{1} = langerman_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = langerman_orig(x);
%
% Purpose:
%
%    Function langerman is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -0.965, for n = 5, 10
%    Global minima: [8.074 8.777 3.467 1.867 6.708]' for n = 5
%                   [8.074 8.777 3.467 1.867 6.708 6.349 4.534...
%                    0.276 7.633 1.567]' for n = 10
%    Local minima 
%    Search domain: 0 <= x <= 10 (Montaz et al. (2005))
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
A   = [0.806 9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020
0.517 9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374
0.100 8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982
0.908 2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426
0.965 8.074 8.777 3.467 1.867 6.708 6.349 4.534 0.276 7.633 1.567];
d = sum((repmat(x',5,1)-A(:,2:11)).^2,2);
f = -sum(A(:,1).*cos(pi.*d).*exp(-d./pi),1);
%
% End of langerman.

