function PB = sphere_k_3D()
%SPHERE_K_3D  DFO Benchmark Suite v2 core: sphere, dimension 3, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K x3:K
%
lb_orig = [-5.1200000000000001 -5.1200000000000001 -5.1200000000000001];
ub_orig = [5.1200000000000001 5.1200000000000001 5.1200000000000001];
lb_C_work = [-5.1200000000000001 -5.1200000000000001 -5.1200000000000001];
ub_C_work = [5.1200000000000001 5.1200000000000001 5.1200000000000001];
T_break = [0 0 0];   % breakpoints of T (NaN for continuous)
T_image = [0 0 0];   % their images
m = 41;
ProblemData = cell(1, 3);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = arrayfun(@(j) sprintf('k3_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'sphere_k_3D';
PB.base_function = 'sphere';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of sphere (F = f o T).';
PB.nvar        = 3;
PB.nobj        = 1;
PB.var_types   = {'K','K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 3);
PB.proportion  = [0 0 3];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [41 41 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 68921;
PB.domain_size = 68921;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'sphere_k_3D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 3);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([38 40 3 2 9 17 8 27 28 41 26 31 19 35 21 7 12 5 13 30 37 39 36 15 10 6 4 18 11 34 24 25 29 20 23 32 22 1 16 33 14])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([21 2 35 22 17 27 31 32 6 13 30 8 1 36 16 9 37 7 15 3 39 38 41 20 14 12 33 40 18 29 4 34 11 5 10 26 23 24 19 25 28])';
g = linspace(lb_orig(3), ub_orig(3), m); D.E{3} = g([6 12 34 3 36 28 31 20 39 14 8 10 2 13 29 11 16 38 37 5 22 40 18 24 19 35 32 4 26 7 15 17 27 25 33 21 23 9 30 41 1])';
D.fn = @(xh) sphere_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [40 28 3];   % protocol starting point: Halton point 1055 (index representation)
PB.x0_alt_idx  = [20 30 31;30 32 12;10 28 39;35 29 19;15 32 9;25 29 27;5 30 17;37 31 35;18 28 25];   % 9 alternative starting points (Halton indices 1055 + 63 k)
PB.core_pidx   = 55;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('sphere_k_3D: n = 3, pattern K, m = 41, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('sphere_k_3D: z must be a finite numeric vector with 3 entries.'); end
xh = zeros(3, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.E{3}(k);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
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
xr = cell(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 3, error('sphere_k_3D: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('sphere_k_3D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('sphere_k_3D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('sphere_k_3D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('sphere_k_3D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = sphere_orig(x)
%
% Purpose:
%
%    Function sphere is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (equal to global minima)
%    Search domain: -5.12 <= x <= 5.12 (Huyer and Neumaier (1999))
%    Cases considered: n = 3 Huyer and Neumaier (1999)
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
f = sum(x.^2,1);
end
%
% End of sphere.

