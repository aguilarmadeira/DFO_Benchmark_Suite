function PB = rastrigin_ck_10D()
%RASTRIGIN_CK_10D  DFO Benchmark Suite v2 core: rastrigin, dimension 10, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K
%
lb_orig = [-5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001];
ub_orig = [5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001];
lb_C_work = [-5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001];
ub_C_work = [5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001];
T_break = [NaN 0 NaN 0 NaN 0 NaN 0 NaN 0];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 NaN 0 NaN 0 NaN 0 NaN 0];   % their images
m = 17;
ProblemData = cell(1, 10);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
ProblemData{10} = arrayfun(@(j) sprintf('k10_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'rastrigin_ck_10D';
PB.base_function = 'rastrigin';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of rastrigin (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C','K','C','K','C','K'};
PB.type_counts = struct('C', 5, 'D', 0, 'K', 5);
PB.proportion  = [5 0 5];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1 1 1];
PB.cardinality = [Inf 17 Inf 17 Inf 17 Inf 17 Inf 17];
PB.n_finite    = 5;
PB.finite_part_size = 1419857;
PB.M_K = 1419857;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'rastrigin_ck_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([14 3 15 13 9 4 2 11 1 6 8 10 5 12 16 17 7])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([10 9 16 3 17 2 12 6 1 8 11 7 5 13 15 14 4])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([11 12 15 16 14 7 8 13 10 9 4 3 5 2 17 1 6])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([17 10 5 15 13 7 16 14 9 1 6 11 4 3 12 8 2])';
g = linspace(lb_orig(10), ub_orig(10), m); D.E{10} = g([11 14 13 12 2 6 8 1 17 15 10 7 4 5 3 16 9])';
D.fn = @(xh) rastrigin_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-4.7949999999999999 14 -3.7568512000000003 10 1.0347708489857261 1 -2.9523631182576837 14 -2.4672113092791985 15];   % protocol starting point: Halton point 1040 (index representation)
PB.x0_alt_idx  = [4.5650000000000004 15 1.4041088000000013 11 -2.1734034560480846 15 4.3821697537146358 2 4.6958231281334752 2;-0.59499999999999975 14 -1.6924671999999998 12 -4.4583621337340347 13 0.91186647669448462 7 2.082590613955781 5;2.0449999999999999 14 3.4684928000000008 10 2.4888354620586028 10 -1.9581681253816408 12 -0.53064190022191138 7;-3.1749999999999998 15 0.29327360000000002 11 0.20387678437265322 8 -4.8282027274577652 1 -3.1438744143996056 10;3.3050000000000006 14 -4.3761663999999998 12 -3.0882262140564167 5 2.5063301445145534 7 4.4635357935399016 13;-1.855 14 2.3576576000000014 12 4.7821870090840806 3 -0.3637044575615711 12 1.8503032793622092 16;0.78500000000000014 15 -2.3117824000000002 11 1.5740127040502703 17 -3.8340077345817223 16 -0.76292923481548502 3;-4.4649999999999999 14 4.4220416000000009 11 -0.71094597363568113 13 3.5005251373905972 5 -3.3761617489931783 5;3.9349999999999996 15 -0.24739840000000068 12 -3.9191202786694896 11 0.63049053531447274 11 3.7868726884194954 8];   % 9 alternative starting points (Halton indices 1040 + 63 k)
PB.core_pidx   = 40;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('rastrigin_ck_10D: n = 10, pattern CK, m = 17, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('rastrigin_ck_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('rastrigin_ck_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('rastrigin_ck_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('rastrigin_ck_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('rastrigin_ck_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('rastrigin_ck_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
k = check_idx(z(10), D.m, D.name, 10); g = D.E{10}(k);
xh(10) = T_map(g, D.lb(10), D.ub(10), D.Tb(10), D.Ti(10));
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
if ~iscell(xr) || numel(xr) ~= 10, error('rastrigin_ck_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('rastrigin_ck_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('rastrigin_ck_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('rastrigin_ck_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('rastrigin_ck_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = rastrigin_orig(x)
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
end
%
% End of rastrigin.

