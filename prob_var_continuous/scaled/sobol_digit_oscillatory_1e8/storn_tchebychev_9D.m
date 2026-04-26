function varargout = storn_tchebychev_9D(varargin)
%STORN_TCHEBYCHEV_9D  storn_tchebychev 9D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-77305.4362523, 23191.6308757]   (range: 100497.067128)
%   x2   ∈ [-1082030.91042, 324609.273125]   (range: 1406640.18354)
%   x3   ∈ [-403028.08654, 120908.425962]   (range: 523936.512502)
%   x4   ∈ [-6781770535.36, 2034531160.61]   (range: 8816301695.96)
%   x5   ∈ [-2803653283.27, 841095984.98]   (range: 3644749268.25)
%   x6   ∈ [-1195097.54733, 358529.264198]   (range: 1553626.81152)
%   x7   ∈ [-519860.826582, 155958.247975]   (range: 675819.074557)
%   x8   ∈ [-8849981664.23, 2654994499.27]   (range: 11504976163.5)
%   x9   ∈ [-85713.6341235, 25714.0902371]   (range: 111427.724361)
%
% Effective contrast ratio (max range / min range): 114480.71563
%
% USAGE:
%   f = storn_tchebychev_9D(x)          % Evaluate function at point x (9D vector)
%   [lb, ub] = storn_tchebychev_9D(n)   % Get bounds for dimension n (must be 9)
%   info = storn_tchebychev_9D()        % Get complete problem information

nloc = 9;
lb_orig = [-128;-128;-128;-128;-128;-128;-128;-128;-128];
ub_orig = [38.40000000000001;38.40000000000001;38.40000000000001;38.40000000000001;38.40000000000001;38.40000000000001;38.40000000000001;38.40000000000001;38.40000000000001];
lb_work = [-77305.43625226879;-1082030.910416686;-403028.0865397486;-6781770535.356745;-2803653283.267669;-1195097.547325445;-519860.8265819909;-8849981664.229654;-85713.63412351672];
ub_work = [23191.63087568064;324609.2731250059;120908.4259619246;2034531160.607024;841095984.9803009;358529.2641976335;155958.2479745973;2654994499.268897;25714.09023705502];
scale_factors = [603.9487207208499;8453.366487630359;3148.656926091786;52982582.30747456;21903541.27552867;9336.699588480038;4061.412707671804;69140481.75179417;669.6377665899744];
contrast_ratio = 114480.7156297488;

% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.

if nargin == 0
    info.name = mfilename;
    info.problem = 'storn_tchebychev';
    info.source_p = 56;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 128;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = false;
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
varargout{1} = storn_tchebychev_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = storn_tchebychev_orig(x);
%
% Purpose:
%
%    Function storn_tchebychev is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima: [128 0 -256 0 160 0 -32 0 1]', for n = 9
%    Local minima
%    Search domain: -128 <= x <= 128, for n = 9
%                   -32768 <= x <= 32768, for n = 17 (Montaz et al. (2005))
%    Cases considered: n = 9, d = 72.661, m = 60
%                      n = 17, d = 10558.145, m = 100 Montaz et al. (2005)
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
if dim == 9
    d = 72.6066;
    m = 60;
else
    if dim == 17
       d = 10558.145;
       m = 100;
    else
        return
    end
end
auxu = (dim-1-[0:dim-1])';
u    = sum(1.2.^auxu.*x,1);
v    = sum((-1.2).^auxu.*x,1);
if u<d
    p1 = (u-d)^2;
else
    p1 = 0;
end
if v<d
    p2 = (v-d)^2;
else
    p2 = 0;
end
w = sum((repmat(2.*[0:m]./m-1,dim,1).^repmat(auxu,1,m+1)).*repmat(x,1,m+1),1);
for i = 1:m
    if w(i) >1
        paux(i,1) = (w(i)-1)^2;
    else
        if w(i)<-1
            paux(i,1) = (w(i)+1)^2;
        else
            paux(i,1) = 0;
        end
    end
end
f = p1+p2+sum(paux,1);
%
% End of storn_tchebychev.

