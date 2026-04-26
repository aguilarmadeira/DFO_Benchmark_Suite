function varargout = rastrigin_10D(varargin)
%RASTRIGIN_10D  rastrigin 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (HALTON OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [-2562560    , 768768      ]   (range: 3331328     )
%   x2   ∈ [-1283.84    , 385.152     ]   (range: 1668.992    )
%   x3   ∈ [-3841280    , 1152384     ]   (range: 4993664     )
%   x4   ∈ [-644.48     , 193.344     ]   (range: 837.824     )
%   x5   ∈ [-3201920    , 960576      ]   (range: 4162496     )
%   x6   ∈ [-1923.2     , 576.96      ]   (range: 2500.16     )
%   x7   ∈ [-4480640    , 1344192     ]   (range: 5824832     )
%   x8   ∈ [-324.8      , 97.44       ]   (range: 422.24      )
%   x9   ∈ [-2882240    , 864672      ]   (range: 3746912     )
%   x10  ∈ [-1603.52    , 481.056     ]   (range: 2084.576    )
%
% Effective contrast ratio (max range / min range): 13795.0738916
%
% Known global minimum (WORK-space):
%   x* = [0;-2.273736754432321e-13;-4.656612873077393e-10;0;-4.656612873077393e-10;0;0;-5.684341886080801e-14;0;0]
%   f* = 0
%
% USAGE:
%   f = rastrigin_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = rastrigin_10D(n)   % Get bounds for dimension n (must be 10)
%   info = rastrigin_10D()        % Get complete problem information

nloc = 10;
lb_orig = [-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12;-5.12];
ub_orig = [1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536;1.536];
lb_work = [-2562560;-1283.84;-3841280;-644.48;-3201920;-1923.2;-4480640;-324.8;-2882240;-1603.52];
ub_work = [768768.0000000003;385.1520000000002;1152384;193.3440000000001;960576.0000000002;576.9600000000003;1344192.000000001;97.44000000000001;864672.0000000003;481.056];
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
    info.problem = 'rastrigin';
    info.source_p = 40;
    info.dimension = nloc;
    info.strategy = 'halton_oscillatory';
    info.kappa = 1000000;
    info.bound_p = 5.12;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = 0;
    info.x_global_min_orig = [0;0;0;0;0;0;0;0;0;0];
    info.x_global_min_work = [0;-2.273736754432321e-13;-4.656612873077393e-10;0;-4.656612873077393e-10;0;0;-5.684341886080801e-14;0;0];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Rastrigin (n=10): x*=0, f*=0. Ref: Ali et al. (2005).';
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
varargout{1} = rastrigin_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = rastrigin_orig(x);
%
% Purpose:
%
%    Function rastrigin is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -5.12 <= x <= 5.12 (Montaz et al. (2005))
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
f   = 10*dim + sum(x.^2-10.*cos(2*pi.*x),1);
%
% End of rastrigin.

