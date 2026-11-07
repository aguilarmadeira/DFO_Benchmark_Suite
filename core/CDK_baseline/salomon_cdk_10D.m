function PB = salomon_cdk_10D()
%SALOMON_CDK_10D  DFO Benchmark Suite v2 core: salomon, dimension 10, pattern CDK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K x9:D x10:C
%
lb_orig = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub_orig = [100 100 100 100 100 100 100 100 100 100];
lb_C_work = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub_C_work = [100 100 100 100 100 100 100 100 100 100];
T_break = [NaN 0 0 NaN 0 0 NaN 0 0 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 0 NaN 0 0 NaN 0 0 NaN];   % their images
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
PB.name        = 'salomon_cdk_10D';
PB.base_function = 'salomon';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CDK adaptation of salomon (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C','D','K','C','K','D','C'};
PB.type_counts = struct('C', 4, 'D', 3, 'K', 3);
PB.proportion  = [4 3 3];
PB.pattern     = 'CDK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1 1];
PB.cardinality = [Inf 17 17 Inf 17 17 Inf 17 17 Inf];
PB.n_finite    = 6;
PB.finite_part_size = 24137569;
PB.M_K = 4913;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'salomon_cdk_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([14 8 5 10 2 15 4 1 7 13 16 3 17 6 12 9 11])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([3 1 9 4 13 17 7 11 6 16 5 15 14 10 8 2 12])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([17 10 1 9 7 2 5 3 6 14 8 4 16 12 13 11 15])';
D.fn = @(xh) salomon_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [56.34765625 16 13 -91.586838817159517 15 5 -22.369224506411584 16 7 94.53852146459468];   % protocol starting point: Halton point 1043 (index representation)
PB.x0_alt_idx  = [-42.08984375 17 5 -82.840483132028325 10 3 -78.424587828210861 4 2 -70.265283529460007;7.12890625 16 16 -74.094127446897133 6 17 53.104009769997958 10 14 -35.306900651933248;-91.30859375 16 9 -93.336109954185758 1 14 -2.951353551801347 15 10 -0.34851777440650267;87.98828125 17 4 -84.589754269054566 14 12 -59.006716873600652 4 6 34.609865103120256;-10.44921875 16 12 -75.843398583923374 9 8 84.245878282108691 9 2 69.568247980647016;38.76953125 16 7 -99.083715118700539 5 6 28.1905149603094 14 14 -95.235557013407686;-59.66796875 16 15 -90.337359433569347 16 4 -39.588845918990437 3 10 -60.277174135880927;62.79296875 15 10 -81.591003748438155 12 1 -95.644209240789735 8 6 -25.318791258354182;-35.64453125 16 3 -72.844648063306948 7 15 47.608385914919609 13 1 9.6395916191725775];   % 9 alternative starting points (Halton indices 1043 + 63 k)
PB.core_pidx   = 43;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('salomon_cdk_10D: n = 10, pattern CDK, m = 17, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('salomon_cdk_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('salomon_cdk_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('salomon_cdk_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
k = check_idx(z(5), D.m, D.name, 5); g = D.lb(5) + (k - 1) * (D.ub(5) - D.lb(5)) / (D.m - 1);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('salomon_cdk_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
k = check_idx(z(9), D.m, D.name, 9); g = D.lb(9) + (k - 1) * (D.ub(9) - D.lb(9)) / (D.m - 1);
xh(9) = T_map(g, D.lb(9), D.ub(9), D.Tb(9), D.Ti(9));
v = z(10); if v < D.lbw(10) || v > D.ubw(10), error('salomon_cdk_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
xh(10) = map_work_to_orig(v, D.lbw(10), D.ubw(10), D.lb(10), D.ub(10));
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
if ~iscell(xr) || numel(xr) ~= 10, error('salomon_cdk_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('salomon_cdk_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('salomon_cdk_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('salomon_cdk_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('salomon_cdk_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x)
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
end
%
% End of salomon.

