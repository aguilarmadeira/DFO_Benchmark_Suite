function PB = neumaier3_cdk_10D()
%NEUMAIER3_CDK_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern CDK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
%
lb_orig = [-100 -100 -100 -72.698042482299044 -66.211597568472087 -71.1692167641189 -100 -100 -100 -100];
ub_orig = [61.813232421875 69.889549357313427 72.739328 100 100 100 62.208259000730351 69.966651575724555 61.652338292101582 65.622042724178925];
lb_C_work = [-100 -100 -100 -7269.8042482299052 -66.211597568472087 -71.1692167641189 -1 -100 -100 -100];
ub_C_work = [61.813232421875 69.889549357313427 72.739328 10000 100 100 0.62208259000730348 69.966651575724555 61.652338292101582 65.622042724178925];
T_break = [NaN 16.799065183152976 18.758287999999993 NaN 27.282426063793466 25.113467665697996 NaN 27.474988681793405 21.239253719076174 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 18 24 NaN 30 30 NaN 24 18 NaN];   % their images
m = 17;
ProblemData = cell(1, 10);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = {linspace(lb_orig(5), ub_orig(5), m)'};
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{9} = {linspace(lb_orig(9), ub_orig(9), m)'};
ProblemData{10} = [lb_C_work(10); ub_C_work(10)];
PB.name        = 'neumaier3_cdk_10D';
PB.base_function = 'neumaier3';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CDK adaptation of neumaier3 (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C','D','K','C','K','D','C'};
PB.type_counts = struct('C', 4, 'D', 3, 'K', 3);
PB.proportion  = [4 3 3];
PB.pattern     = 'CDK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01 1];
PB.cardinality = [Inf 17 17 Inf 17 17 Inf 17 17 Inf];
PB.n_finite    = 6;
PB.finite_part_size = 24137569;
PB.M_K = 4913;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'neumaier3_cdk_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.cw = [10 NaN NaN 2799.9999999999991 NaN NaN 0.2799999999999998 NaN NaN 10];   % anchors of the continuous maps (NaN for finite)
D.cx = [10 NaN NaN 28 NaN NaN 28 NaN NaN 10];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([6 5 10 12 14 3 16 11 17 1 2 7 4 15 8 13 9])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([3 6 1 10 2 15 14 16 9 13 8 11 4 17 5 7 12])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([16 15 4 6 12 2 9 10 8 5 1 7 3 13 11 17 14])';
D.fn = @(xh) neumaier3_orig(xh);
PB.f_star      = -210;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [31.552261710166931 2 2 7554.4633717625302 3 11 0.48275451083305509 9 2 15.41049719535124];   % protocol starting point: Halton point 1035 (index representation)
PB.x0_alt_idx  = [-48.09018862247467 2 14 8309.7026246005735 14 9 0.02812236622893205 14 13 44.359891107900921;-8.2689634561538696 1 6 9064.941877438614 11 7 -0.52159623999173976 3 9 -92.11582305126197;-87.911413788795471 2 17 9820.1811302766546 6 4 -0.97622838459586281 8 5 -63.166429138712267;41.981630206108093 3 9 8158.6547740329652 2 2 0.18560931828134075 13 17 -34.217035226162565;-37.660820126533508 2 2 8913.8940268710066 13 16 -0.2690228263227824 2 12 -5.2676413136128417;2.1604050397872925 2 12 9323.8810498402272 9 13 -0.72365497092690556 8 8 23.681752598936853;-77.482045292854309 2 5 7662.3546935965342 4 11 0.34309627033374945 12 4 52.631146511486548;21.596955418586731 1 15 8417.5939464345756 16 7 -0.11153587427037359 1 16 -83.844567647676342;-58.045494914054871 2 8 9172.8331992726162 11 5 -0.56616801887449697 6 11 -54.89517373512664];   % 9 alternative starting points (Halton indices 1035 + 63 k)
PB.core_pidx   = 35;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('neumaier3_cdk_10D: n = 10, pattern CDK, m = 17, f* = -210 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('neumaier3_cdk_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('neumaier3_cdk_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('neumaier3_cdk_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
k = check_idx(z(5), D.m, D.name, 5); g = D.lb(5) + (k - 1) * (D.ub(5) - D.lb(5)) / (D.m - 1);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('neumaier3_cdk_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
k = check_idx(z(9), D.m, D.name, 9); g = D.lb(9) + (k - 1) * (D.ub(9) - D.lb(9)) / (D.m - 1);
xh(9) = T_map(g, D.lb(9), D.ub(9), D.Tb(9), D.Ti(9));
v = z(10); tl = 1e-12 * (D.ubw(10) - D.lbw(10)); if v < D.lbw(10) - tl || v > D.ubw(10) + tl, error('neumaier3_cdk_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
v = min(max(v, D.lbw(10)), D.ubw(10));
if abs(v - D.cw(10)) <= 1e-12 * (D.ubw(10) - D.lbw(10)), xh(10) = D.cx(10); else, xh(10) = map_work_to_orig(v, D.lbw(10), D.ubw(10), D.lb(10), D.ub(10)); end
end

function x = T_map(g, a, b, zs, xs)
% piecewise-affine T with T(zs) = xs; fixes the bounds when zs is interior
if abs(g - zs) <= 1e-12 * max(1, abs(b - a)), x = xs; return; end
if abs(zs - a) <= 1e-14 * max(1, abs(a)), x = xs + (g - a) * (b - xs) / (b - a);
elseif abs(zs - b) <= 1e-14 * max(1, abs(b)), x = a + (g - a) * (xs - a) / (b - a);
elseif g <= zs, x = a + (g - a) * (xs - a) / (zs - a);
else, x = xs + (g - zs) * (b - xs) / (b - zs);
end
x = min(max(x, min(a, xs)), max(b, xs));
end

function k = check_idx(v, m, wname, i)
if ~isscalar(v) || ~isfinite(v) || v ~= floor(v) || v < 1 || v > m, error('%s: z%d = %g is not an index in 1..%d.', wname, i, v, m); end
k = double(v);
end

function xo = map_work_to_orig(xw, lbw, ubw, lbo, ubo)
if abs(ubw - lbw) <= 1e-14, t = 0; else, t = (xw - lbw) / (ubw - lbw); end
t = max(0, min(1, t)); xo = lbo + t * (ubo - lbo);
end

function xr = idx2cell(z, PB)
xr = cell(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 10, error('neumaier3_cdk_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('neumaier3_cdk_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('neumaier3_cdk_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('neumaier3_cdk_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('neumaier3_cdk_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = neumaier3_orig(x)
%
% Purpose:
%
%    Function neumaier3 is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -n*(n+4)*(n-1)/6
%    Global minima = [1:n]'.*(n+1-[1:n]')
%    Local minima (several)
%    Search domain: -n^2 <= x <= n^2   (Montaz et al. (2005))
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
dim  = size(x,1);
xaux = x(1:dim-1,1);
f    = sum((x-1).^2,1)-sum(x(2:dim,1).*xaux,1);
end
%
% End of neumaier3.

