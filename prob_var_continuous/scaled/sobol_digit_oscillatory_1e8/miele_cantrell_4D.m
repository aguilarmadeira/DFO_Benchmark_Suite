function varargout = miele_cantrell_4D(varargin)
%MIELE_CANTRELL_4D  miele_cantrell 4D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-903108011.996, 270932403.599]   (range: 1174040415.6)
%   x2   ∈ [-362493620.677, 108748086.203]   (range: 471241706.88)
%   x3   ∈ [-60264.1098173, 18079.2329452]   (range: 78343.3427625)
%   x4   ∈ [-13188.1906607, 3956.4571982]   (range: 17144.6478589)
%
% Effective contrast ratio (max range / min range): 68478.5377489
%
% Known global minimum (WORK-space):
%   x* = [1.192092895507812e-07;36249362.06766272;6026.410981731737;1318.819066066717]
%   f* = 0
%
% USAGE:
%   f = miele_cantrell_4D(x)          % Evaluate function at point x (4D vector)
%   [lb, ub] = miele_cantrell_4D(n)   % Get bounds for dimension n (must be 4)
%   info = miele_cantrell_4D()        % Get complete problem information

nloc = 4;
lb_orig = [-10;-10;-10;-10];
ub_orig = [3;3;3;3];
lb_work = [-903108011.9964733;-362493620.6766272;-60264.10981731743;-13188.19066066716];
ub_work = [270932403.598942;108748086.2029881;18079.23294519523;3956.457198200149];
scale_factors = [90310801.19964734;36249362.06766272;6026.410981731742;1318.819066066716];
contrast_ratio = 68478.53774892175;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'miele_cantrell';
    info.source_p = 32;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 10;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;1;1;1];
    info.x_global_min_work = [1.192092895507812e-07;36249362.06766272;6026.410981731737;1318.819066066717];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Miele-Cantrell (4D): x*=(0,1,1,1), f*=0. Ref: Brachetti et al. (1997).';
    info.mapping = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1} = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 4D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 4 components.');
end
range = ub_work - lb_work;
range(range == 0) = 1;
t = (x - lb_work)./range;
t = max(0, min(1, t));
x_orig = lb_orig + t.*(ub_orig - lb_orig);
varargout{1} = miele_cantrell_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = miele_cantrell_orig(x);
%
% Purpose:
%
%    Function miele_cantrell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [0 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = (exp(x(1))-x(2))^4 + 100*(x(2)-x(3))^6 + tan(x(3)-x(4))^4 + x(1)^8;
%
% End of miele_cantrell.

