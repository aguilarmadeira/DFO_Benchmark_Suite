function varargout = griewank_10D(varargin)
%GRIEWANK_10D  griewank 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-20059573617.6, 6017872085.27]   (range: 26077445702.8)
%   x2   ∈ [-1765210.28306, 529563.084917]   (range: 2294773.36797)
%   x3   ∈ [-36893133697.3, 11067940109.2]   (range: 47961073806.5)
%   x4   ∈ [-449020.495183, 134706.148555]   (range: 583726.643738)
%   x5   ∈ [-27770233387.4, 8331070016.21]   (range: 36101303403.6)
%   x6   ∈ [-10171064748.9, 3051319424.66]   (range: 13222384173.5)
%   x7   ∈ [-3468550.1605, 1040565.04815]   (range: 4509115.20865)
%   x8   ∈ [-708338.790614, 212501.637184]   (range: 920840.427798)
%   x9   ∈ [-2481786.22832, 744535.868495]   (range: 3226322.09681)
%   x10  ∈ [-17263484803.7, 5179045441.11]   (range: 22442530244.8)
%
% Effective contrast ratio (max range / min range): 82163.5851662
%
% Known global minimum (WORK-space):
%   x* = [3.814697265625e-06;0;0;5.820766091346741e-11;3.814697265625e-06;1.9073486328125e-06;4.656612873077393e-10;0;0;3.814697265625e-06]
%   f* = 0
%
% USAGE:
%   f = griewank_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = griewank_10D(n)   % Get bounds for dimension n (must be 10)
%   info = griewank_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-400;-400;-400;-400;-400;-400;-400;-400;-400;-400];
ub_orig = [120;120;120;120;120;120;120;120;120;120];
lb_work = [-20059573617.55868;-1765210.283056742;-36893133697.32005;-449020.4951829637;-27770233387.3513;-10171064748.87625;-3468550.160498296;-708338.7906142104;-2481786.228316848;-17263484803.70806];
ub_work = [6017872085.267604;529563.0849170225;11067940109.19601;134706.1485548891;8331070016.205389;3051319424.662876;1040565.048149489;212501.6371842631;744535.8684950545;5179045441.112418];
scale_factors = [50148934.0438967;4413.025707641855;92232834.24330011;1122.551237957409;69425583.46837825;25427661.87219063;8671.375401245739;1770.846976535526;6204.46557079212;43158712.00927015];
contrast_ratio = 82163.58516616728;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'griewank';
    info.source_p = 22;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 400;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [3.814697265625e-06;0;0;5.820766091346741e-11;3.814697265625e-06;1.9073486328125e-06;4.656612873077393e-10;0;0;3.814697265625e-06];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Griewank (n=10): x*=0, f*=0. Ref: Huyer & Neumaier (1999).';
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

