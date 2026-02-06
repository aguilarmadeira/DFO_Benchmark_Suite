function varargout = paviani_10D(varargin)
%PAVIANI_10D  paviani 10D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (SOBOL DIGIT OSCILLATORY HETEROGENEITY):
%
%   x1   ∈ [17258.7455127, 86241.9772021]   (range: 68983.2316894)
%   x2   ∈ [11783613.9685, 58882736.667]   (range: 47099122.6986)
%   x3   ∈ [11562.024472, 57775.4536207]   (range: 46213.4291488)
%   x4   ∈ [50167834.0571, 250688741.997]   (range: 200520907.94)
%   x5   ∈ [19523.5615902, 97559.2665368]   (range: 78035.7049467)
%   x6   ∈ [32869463.7458, 164248759.617]   (range: 131379295.871)
%   x7   ∈ [14295.5442825, 71434.8562122]   (range: 57139.3119297)
%   x8   ∈ [7906.53491111, 39508.9668047]   (range: 31602.4318936)
%   x9   ∈ [15031.1419259, 75110.6387394]   (range: 60079.4968134)
%   x10  ∈ [15293104.9689, 76419668.4579]   (range: 61126563.489)
%
% Effective contrast ratio (max range / min range): 6345.11004139
%
% Known global minimum (WORK-space):
%   x* = [80652.93817551536;55066753.73270464;54031.22980372891;234442486.8904013;91236.79381797538;153604375.5556635;66805.41458546891;36948.52971203888;70242.98258344817;71467178.69285187]
%   f* = -45.778
%
% USAGE:
%   f = paviani_10D(x)          % Evaluate function at point x (10D vector)
%   [lb, ub] = paviani_10D(n)   % Get bounds for dimension n (must be 10)
%   info = paviani_10D()        % Get complete problem information

nloc = 10;
lb_orig = [2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001;2.001];
ub_orig = [9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001;9.999000000000001];
lb_work = [17258.74551269449;11783613.96846775;11562.0244719561;50167834.05707334;19523.56159017952;32869463.74579003;14295.54428248564;7906.534911110017;15031.14192594158;15293104.96892274];
ub_work = [86241.97720211509;58882736.66702103;57775.45362073419;250688741.9973396;97559.2665368341;164248759.6172687;71434.85621218089;39508.96680469219;75110.63873937528;76419668.45790032];
scale_factors = [8625.060226234131;5888862.552957399;5778.123174390857;25071381.33786775;9756.9023439178;16426518.61358823;7144.200041222211;3951.291809650184;7511.815055443072;7642731.118901921];
contrast_ratio = 6345.110041388558;

% Reference:
%   J. F. A. Madeira,
%   "Global and Local Optimization using Direct Search - A Scale-Invariant Approach (GLODS-SI)",
%   2026.

if nargin == 0
    info.name = mfilename;
    info.problem = 'paviani';
    info.source_p = 36;
    info.dimension = nloc;
    info.strategy = 'sobol_digit_oscillatory';
    info.kappa = 100000000;
    info.bound_p = 0;  % Original bound(p) from test setup
    info.lb_orig = lb_orig; info.ub_orig = ub_orig;
    info.lb_work = lb_work; info.ub_work = ub_work;
    info.scale_factors = scale_factors;
    info.contrast_ratio = contrast_ratio;
    info.global_min_known = true;
    info.f_global_min = -45.778;
    info.x_global_min_orig = [9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001;9.351000000000001];
    info.x_global_min_work = [80652.93817551536;55066753.73270464;54031.22980372891;234442486.8904013;91236.79381797538;153604375.5556635;66805.41458546891;36948.52971203888;70242.98258344817;71467178.69285187];
    info.global_min_note = 'Mapped x*_orig -> x*_work via affine inverse using t=(x*_orig-lb_orig)./(ub_orig-lb_orig). Original note: Paviani (10D): x*=9.351, f*=-45.778. Ref: Ali et al. (2005).';
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
varargout{1} = paviani_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original problem function (verbatim; only the function name is renamed)
% -------------------------------------------------------------------------
function [f] = paviani_orig(x);
%
% Purpose:
%
%    Function paviani is the function described in
%    Montaz et al. (2005).
%
%    dim = 10
%    Minimum global value = -45.778
%    Global minima = 9.351*ones(10,1)
%    Local minima
%    Search domain: 2.001 <= x <= 9.999
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
f = sum((log(x-2)).^2+(log(10-x)).^2,1)-prod(x,1)^0.2;
%
% End of paviani.

