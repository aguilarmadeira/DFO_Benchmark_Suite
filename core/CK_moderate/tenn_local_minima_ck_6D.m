function PB = tenn_local_minima_ck_6D()
%TENN_LOCAL_MINIMA_CK_6D  DFO Benchmark Suite v2 core: tenn_local_minima, dimension 6, pattern CK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 7.852723858933194e-32 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C x6:K
%
lb_orig = [-10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10];
lb_C_work = [-10 -10 -1000 -10 -0.10000000000000001 -10];
ub_C_work = [10 10 1000 10 0.10000000000000001 10];
T_break = [NaN 1 NaN 1 NaN 1];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 NaN 1 NaN 1];   % their images
m = 41;
ProblemData = cell(1, 6);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'tenn_local_minima_ck_6D';
PB.base_function = 'tenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of tenn_local_minima (F = f o T).';
PB.nvar        = 6;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C','K'};
PB.type_counts = struct('C', 3, 'D', 0, 'K', 3);
PB.proportion  = [3 0 3];
PB.pattern     = 'CK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01];
PB.cardinality = [Inf 41 Inf 41 Inf 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 68921;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'tenn_local_minima_ck_6D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 6);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([15 6 12 19 30 11 10 29 24 28 17 38 36 16 21 37 27 14 23 1 5 39 9 25 4 31 32 41 35 13 34 18 7 20 3 22 33 2 8 40 26])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([26 8 24 23 27 20 16 28 6 2 19 4 41 22 40 32 30 29 14 7 31 15 17 37 35 33 3 12 18 1 36 38 5 25 13 11 34 10 21 9 39])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([28 8 34 23 10 9 32 18 7 21 11 6 38 15 20 14 33 29 30 12 17 39 1 40 22 41 16 26 5 19 2 3 4 37 24 36 31 13 35 25 27])';
D.fn = @(xh) tenn_local_minima_orig(xh);
PB.f_star      = 7.8527238589331936e-32;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [5.322265625 10 722.24000000000024 16 -0.031029301277235175 20];   % protocol starting point: Halton point 1059 (index representation)
PB.x0_alt_idx  = [-4.521484375 13 194.24000000000024 17 -0.093688955672426757 15;0.400390625 14 -794.55999999999995 14 0.060030052592036059 7;-9.443359375 11 597.44000000000028 15 -0.0026296018031555096 2;9.423828125 12 -406.72000000000003 17 -0.047257700976709241 37;-0.419921875 14 985.2800000000002 12 0.086790519773239566 29;4.501953125 11 -3.5199999999999818 14 0.042162420599685807 24;-5.341796875 12 -992.32000000000005 16 -0.020497233795505776 19;6.904296875 13 399.68000000000006 18 -0.065125332969059493 11;-2.939453125 10 -589.12 14 0.070562120073765472 6];   % 9 alternative starting points (Halton indices 1059 + 63 k)
PB.core_pidx   = 59;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('tenn_local_minima_ck_6D: n = 6, pattern CK, m = 41, f* = 7.852723859e-32 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 6 || ~all(isfinite(z)), error('tenn_local_minima_ck_6D: z must be a finite numeric vector with 6 entries.'); end
xh = zeros(6, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('tenn_local_minima_ck_6D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('tenn_local_minima_ck_6D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('tenn_local_minima_ck_6D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
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
xr = cell(1, 6);
for i = 1:6
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 6, error('tenn_local_minima_ck_6D: xr must be a cell array with 6 entries.'); end
z = zeros(1, 6);
for i = 1:6
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_ck_6D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('tenn_local_minima_ck_6D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('tenn_local_minima_ck_6D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('tenn_local_minima_ck_6D: x%d unknown label %s.', i, lab); end
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

