function varargout = storn_tchebychev_9D(varargin)
%STORN_TCHEBYCHEV_9D  storn_tchebychev 9D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-64064000   , 64064000    ]   (range: 128128000   )
%   x2   ∈ [-32096      , 32096       ]   (range: 64192       )
%   x3   ∈ [-96032000   , 96032000    ]   (range: 192064000   )
%   x4   ∈ [-16112      , 16112       ]   (range: 32224       )
%   x5   ∈ [-80048000   , 80048000    ]   (range: 160096000   )
%   x6   ∈ [-48080      , 48080       ]   (range: 96160       )
%   x7   ∈ [-112016000  , 112016000   ]   (range: 224032000   )
%   x8   ∈ [-8120       , 8120        ]   (range: 16240       )
%   x9   ∈ [-72056000   , 72056000    ]   (range: 144112000   )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% USAGE:
%   f = storn_tchebychev_9D(x)          % Evaluate function at point x (9D vector)
%   [lb, ub] = storn_tchebychev_9D(n)   % Get bounds for dimension n (must be 9)
%   info = storn_tchebychev_9D()        % Get complete problem information

nloc = 9;
lb_orig = [-128;-128;-128;-128;-128;-128;-128;-128;-128];
ub_orig = [128;128;128;128;128;128;128;128;128];
lb_work = [-64064000;-32096;-96032000;-16112;-80048000;-48080;-112016000;-8120;-72056000];
ub_work = [64064000;32096;96032000;16112;80048000;48080;112016000;8120;72056000];
scale_factors = [500500;250.75;750250;125.875;625375;375.625;875125;63.4375;562937.5];
contrast_ratio = 13795.07389162562;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'storn_tchebychev';
    info.source_p = 56;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
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

