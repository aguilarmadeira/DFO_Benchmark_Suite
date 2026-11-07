function PB = salomon_k_10D()
%SALOMON_K_10D  DFO Benchmark Suite v2 core: salomon, dimension 10, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 9 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K x3:K x4:K x5:K x6:K x7:K x8:K x9:K x10:K
%
lb_orig = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub_orig = [100 100 100 100 100 100 100 100 100 100];
lb_C_work = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub_C_work = [100 100 100 100 100 100 100 100 100 100];
T_break = [0 0 0 0 0 0 0 0 0 0];   % breakpoints of T (NaN for continuous)
T_image = [0 0 0 0 0 0 0 0 0 0];   % their images
m = 9;
ProblemData = cell(1, 10);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = arrayfun(@(j) sprintf('k3_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = arrayfun(@(j) sprintf('k5_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = arrayfun(@(j) sprintf('k7_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{9} = arrayfun(@(j) sprintf('k9_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{10} = arrayfun(@(j) sprintf('k10_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'salomon_k_10D';
PB.base_function = 'salomon';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of salomon (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'K','K','K','K','K','K','K','K','K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 10);
PB.proportion  = [0 0 10];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [9 9 9 9 9 9 9 9 9 9];
PB.n_finite    = 10;
PB.finite_part_size = 3486784401;
PB.M_K = 3486784401;
PB.domain_size = 3486784401;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'salomon_k_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([7 3 1 2 8 9 4 6 5])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([4 2 7 1 8 3 9 5 6])';
g = linspace(lb_orig(3), ub_orig(3), m); D.E{3} = g([8 5 9 1 4 7 6 3 2])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([7 6 1 5 4 3 2 8 9])';
g = linspace(lb_orig(5), ub_orig(5), m); D.E{5} = g([6 7 5 9 8 1 2 4 3])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([5 8 2 3 9 7 4 1 6])';
g = linspace(lb_orig(7), ub_orig(7), m); D.E{7} = g([3 5 9 8 4 2 6 7 1])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([3 1 5 2 4 7 9 6 8])';
g = linspace(lb_orig(9), ub_orig(9), m); D.E{9} = g([3 6 5 8 9 7 2 1 4])';
g = linspace(lb_orig(10), ub_orig(10), m); D.E{10} = g([1 8 3 9 4 6 2 7 5])';
D.fn = @(xh) salomon_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [7 9 7 1 8 3 4 9 4 9];   % protocol starting point: Halton point 1043 (index representation)
PB.x0_alt_idx  = [3 9 3 2 5 2 2 3 2 2;5 8 8 2 4 9 7 5 8 4;1 8 5 1 1 7 5 8 6 5;9 9 2 2 7 6 3 3 4 6;5 8 6 2 5 5 8 5 2 8;7 9 4 1 3 4 6 8 7 1;3 9 8 1 8 3 3 2 5 3;8 8 6 2 7 1 1 4 3 4;4 9 2 2 4 8 7 7 1 5];   % 9 alternative starting points (Halton indices 1043 + 63 k)
PB.core_pidx   = 43;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('salomon_k_10D: n = 10, pattern K, m = 9, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('salomon_k_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.E{3}(k);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
k = check_idx(z(5), D.m, D.name, 5); g = D.E{5}(k);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
k = check_idx(z(7), D.m, D.name, 7); g = D.E{7}(k);
xh(7) = T_map(g, D.lb(7), D.ub(7), D.Tb(7), D.Ti(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
k = check_idx(z(9), D.m, D.name, 9); g = D.E{9}(k);
xh(9) = T_map(g, D.lb(9), D.ub(9), D.Tb(9), D.Ti(9));
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
if ~iscell(xr) || numel(xr) ~= 10, error('salomon_k_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('salomon_k_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('salomon_k_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('salomon_k_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('salomon_k_10D: x%d unknown label %s.', i, lab); end
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

