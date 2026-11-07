function PB = tenn_local_minima_ck_8D()
%TENN_LOCAL_MINIMA_CK_8D  DFO Benchmark Suite v2 core: tenn_local_minima, dimension 8, pattern CK, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 25 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 5.889542894199896e-32 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10 10 10];
lb_C_work = [-10 -10 -10 -10 -1000000000 -10 -1000000000 -10];
ub_C_work = [10 10 10 10 1000000000 10 1000000000 10];
T_break = [NaN 0.83333333333333337 NaN 0.83333333333333337 NaN 0.83333333333333337 NaN 0.83333333333333337];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 NaN 1 NaN 1 NaN 1];   % their images
m = 25;
ProblemData = cell(1, 8);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'tenn_local_minima_ck_8D';
PB.base_function = 'tenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of tenn_local_minima (F = f o T).';
PB.nvar        = 8;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C','K','C','K'};
PB.type_counts = struct('C', 4, 'D', 0, 'K', 4);
PB.proportion  = [4 0 4];
PB.pattern     = 'CK';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 1 100000000 100000000];
PB.cardinality = [Inf 25 Inf 25 Inf 25 Inf 25];
PB.n_finite    = 4;
PB.finite_part_size = 390625;
PB.M_K = 390625;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'tenn_local_minima_ck_8D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 8);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([4 11 19 14 3 5 10 25 7 23 17 2 6 22 18 21 24 12 16 1 15 20 13 8 9])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([18 6 4 23 5 8 12 1 16 14 19 9 13 24 17 11 7 20 25 21 15 3 22 10 2])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([17 12 24 13 1 2 11 14 25 5 19 15 21 22 3 23 6 18 8 9 16 4 7 20 10])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([10 17 11 13 16 4 25 7 8 23 2 22 24 12 19 18 5 20 21 3 14 9 15 1 6])';
D.fn = @(xh) tenn_local_minima_orig(xh);
PB.f_star      = 5.889542894199896e-32;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-7.177734375 14 -7.9775999999999998 13 -128474830.95416975 14 -216771829.83920217 21];   % protocol starting point: Halton point 1060 (index representation)
PB.x0_alt_idx  = [5.478515625 16 5.9424000000000028 14 -755071374.90608561 11 -777325463.05719519 4;-4.599609375 17 -3.9455999999999989 12 782118707.73854256 7 537960512.92489314 11;0.556640625 15 9.9744000000000028 13 155522163.78662658 4 -22593120.293099999 19;-9.638671875 15 -0.067199999999999704 14 -290758827.94891071 25 -583146753.51109302 3;9.580078125 17 -9.9551999999999996 11 -933747694.82958817 20 849379198.04600048 11;-0.498046875 15 3.9648000000000021 12 603442387.81503987 17 288825564.82800746 19;4.658203125 16 -5.9231999999999996 13 -23154156.136875868 14 -388968043.96499085 1;-5.595703125 16 7.9968000000000004 14 -469435147.87241316 9 -949521677.18298388 9;7.060546875 15 -1.8911999999999995 12 887439382.55583644 6 483004274.37410951 17];   % 9 alternative starting points (Halton indices 1060 + 63 k)
PB.core_pidx   = 60;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('tenn_local_minima_ck_8D: n = 8, pattern CK, m = 25, f* = 5.889542894e-32 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 8 || ~all(isfinite(z)), error('tenn_local_minima_ck_8D: z must be a finite numeric vector with 8 entries.'); end
xh = zeros(8, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('tenn_local_minima_ck_8D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('tenn_local_minima_ck_8D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('tenn_local_minima_ck_8D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('tenn_local_minima_ck_8D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
end

function x = T_map(g, a, b, zs, xs)
% piecewise-affine T with T(zs) = xs; fixes the bounds when zs is interior
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
xr = cell(1, 8);
for i = 1:8
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 8, error('tenn_local_minima_ck_8D: xr must be a cell array with 8 entries.'); end
z = zeros(1, 8);
for i = 1:8
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_ck_8D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('tenn_local_minima_ck_8D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('tenn_local_minima_ck_8D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('tenn_local_minima_ck_8D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = tenn_local_minima_orig(x)
%
% Purpose:
%
%    Function tenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (10^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
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
f   = pi/dim*sum(((x(1:dim-1)-1).^2).*(1+10*sin(pi.*x(2:dim)).^2),1);
f   = f + pi/dim*(10*sin(pi*x(1))^2+(x(dim)-1)^2);
end
%
% End of tenn_local_minima.

