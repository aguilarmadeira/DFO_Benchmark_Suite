function PB = kowalik_k_4D()
%KOWALIK_K_4D  DFO Benchmark Suite v2.1.0 core: kowalik, dimension 4, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0.0003074859878056062 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K x3:K x4:K
%
lb_orig = [0.059979711914062504 0.081743057460752944 0 0.056734110787172012];
ub_orig = [0.41999999999999998 0.41999999999999998 0.34746458879999997 0.41999999999999998];
lb_C_work = [0.059979711914062504 0.081743057460752944 0 0.056734110787172012];
ub_C_work = [0.41999999999999998 0.41999999999999998 0.34746458879999997 0.41999999999999998];
T_break = [0.19498731994628907 0.19449537164050196 0.13029922079999998 0.13241450437317784];   % breakpoints of T (NaN for continuous)
T_image = [0.1928334529007549 0.1908362411132952 0.1231173000994744 0.1357659904379078];   % their images
m = 25;
ProblemData = cell(1, 4);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = arrayfun(@(j) sprintf('k3_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'kowalik_k_4D';
PB.base_function = 'kowalik';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of kowalik (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'K','K','K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 4);
PB.proportion  = [0 0 4];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [25 25 25 25];
PB.n_finite    = 4;
PB.finite_part_size = 390625;
PB.M_K = 390625;
PB.domain_size = 390625;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'kowalik_k_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [NaN NaN NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([18 2 1 15 19 16 3 6 23 14 5 25 7 22 4 8 17 13 11 9 12 21 10 20 24])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([16 9 19 11 22 18 4 7 2 1 23 21 24 3 25 14 8 20 5 6 12 13 17 15 10])';
g = linspace(lb_orig(3), ub_orig(3), m); D.E{3} = g([11 9 2 14 6 8 5 13 19 12 21 20 22 7 15 4 18 25 17 3 10 24 23 1 16])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([3 13 24 1 7 17 12 22 16 15 23 6 4 11 5 19 8 9 21 20 10 14 25 18 2])';
D.fn = @(xh) kowalik_orig(xh);
PB.f_star      = 0.00030748598780560617;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [4 18 16 25];   % protocol starting point: Halton point 1028 (index representation)
PB.x0_alt_idx  = [19 18 9 22;7 19 21 23;13 17 14 24;1 19 2 22;25 19 19 23;13 18 7 24;19 18 24 25;7 19 12 22;22 18 6 23];   % 9 alternative starting points (Halton indices 1028 + 63 k)
PB.core_pidx   = 28;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('kowalik_k_4D: n = 4, pattern K, m = 25, f* = 0.0003074859878 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('kowalik_k_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.E{3}(k);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
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
xr = cell(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 4, error('kowalik_k_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('kowalik_k_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('kowalik_k_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('kowalik_k_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('kowalik_k_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = kowalik_orig(x)
%
% Purpose:
%
%    Function kowalik is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 3.0748*10^-4
%    Global minima = [0.192 0.190 0.123 0.135]'
%    Local minima 
%    Search domain: 0 <= x <= 0.42 (Montaz et al. (2005))
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
a = [0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b = [0.25 0.5 1 2 4 6 8 10 12 14 16];
f = sum((a-x(1).*(1+x(2).*b)./(1+x(3).*b+x(4).*b.^2)).^2,2);
end
%
% End of kowalik.

