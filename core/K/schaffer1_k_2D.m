function PB = schaffer1_k_2D()
%SCHAFFER1_K_2D  DFO Benchmark Suite v2 core: schaffer1, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-100 -100];
ub_orig = [100 100];
lb_C_work = [-100 -100];
ub_C_work = [100 100];
T_break = [0 0];   % breakpoints of T (NaN for continuous)
T_image = [0 0];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'schaffer1_k_2D';
PB.base_function = 'schaffer1';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of schaffer1 (F = f o T).';
PB.nvar        = 2;
PB.nobj        = 1;
PB.var_types   = {'K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 2);
PB.proportion  = [0 0 2];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [151 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = 22801;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'schaffer1_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([150 42 55 129 71 80 52 35 25 119 60 58 133 118 1 31 139 5 56 114 120 98 36 43 32 51 87 123 89 82 135 100 7 41 131 76 27 77 2 6 132 19 22 29 78 137 115 146 81 117 105 8 40 75 13 26 128 116 4 79 102 70 108 20 44 63 10 142 37 85 73 99 97 96 65 18 38 69 126 68 127 141 62 45 148 24 57 74 21 130 112 3 47 23 48 28 107 95 50 83 84 54 124 109 103 67 11 145 140 30 94 72 64 106 125 53 33 147 88 101 9 15 134 91 136 111 86 34 104 61 16 93 143 17 39 12 144 110 46 113 49 138 66 59 121 151 90 92 122 14 149])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([117 81 128 63 87 43 129 35 44 28 116 85 113 29 102 68 124 46 33 75 110 119 57 23 104 127 12 107 115 151 98 150 8 53 139 18 56 1 3 138 47 137 90 74 38 59 78 16 112 148 82 9 134 66 14 86 131 54 62 120 118 109 10 22 45 83 11 31 39 5 106 52 13 123 67 20 136 80 25 71 126 93 65 48 147 94 60 77 24 143 76 6 95 122 51 32 108 2 55 133 89 58 96 69 4 73 17 140 145 37 146 130 72 132 149 61 42 141 97 111 19 27 30 103 26 125 100 144 121 79 36 114 88 41 64 40 91 7 84 15 99 50 135 142 70 101 92 34 21 105 49])';
D.fn = @(xh) schaffer1_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [25 16];   % protocol starting point: Halton point 1044 (index representation)
PB.x0_alt_idx  = [119 6;44 10;83 13;6 2;143 11;68 16;106 3;29 7;124 17];   % 9 alternative starting points (Halton indices 1044 + 63 k)
PB.core_pidx   = 44;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('schaffer1_k_2D: n = 2, pattern K, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('schaffer1_k_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
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
xr = cell(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 2, error('schaffer1_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('schaffer1_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('schaffer1_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('schaffer1_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('schaffer1_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = schaffer1_orig(x)
%
% Purpose:
%
%    Function schaffer1 is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
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
f = 0.5 + (sin(norm(x,2))^2-0.5)/((1+0.001*(x(1)^2+x(2)^2))^2);
end
%
% End of schaffer1.

