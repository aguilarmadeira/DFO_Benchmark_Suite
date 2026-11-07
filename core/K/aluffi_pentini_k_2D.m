function PB = aluffi_pentini_k_2D()
%ALUFFI_PENTINI_K_2D  DFO Benchmark Suite v2 core: aluffi_pentini, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -0.3523860738000364 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-10 -10];
ub_orig = [10 10];
lb_C_work = [-10 -10];
ub_C_work = [10 10];
T_break = [-1.0666666666666667 0];   % breakpoints of T (NaN for continuous)
T_image = [-1.0466805324973201 0];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'aluffi_pentini_k_2D';
PB.base_function = 'aluffi_pentini';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of aluffi_pentini (F = f o T).';
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
D = struct('name', 'aluffi_pentini_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([106 48 86 11 32 128 95 14 22 43 91 59 7 77 142 139 144 141 135 21 114 3 44 47 124 12 79 88 103 23 131 100 99 125 4 5 112 121 146 56 24 109 108 67 80 83 90 122 87 115 27 33 51 54 13 123 64 151 132 52 129 16 19 97 39 50 58 2 6 145 30 118 72 111 150 134 57 82 15 40 61 102 1 38 62 37 116 65 127 69 137 143 74 113 34 76 25 73 60 78 35 75 96 36 71 17 46 105 147 84 9 138 20 107 18 55 110 53 119 148 120 8 94 81 126 68 85 63 93 117 42 41 92 45 136 31 104 149 133 98 130 70 49 140 26 66 28 89 10 29 101])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([50 142 72 131 125 126 111 22 87 35 71 6 134 85 135 49 145 95 96 67 46 53 146 47 105 12 15 141 138 123 51 147 33 100 88 130 10 132 18 56 25 97 75 36 86 2 7 119 13 19 70 109 66 90 14 21 115 37 42 54 74 128 122 113 84 92 60 151 8 28 32 55 133 11 139 48 137 98 112 43 143 93 108 104 65 129 58 68 89 73 103 107 78 24 38 102 29 77 116 83 40 144 101 41 76 23 149 3 34 91 16 17 27 117 26 64 1 63 121 127 9 20 106 94 81 82 99 79 110 5 61 62 114 69 57 59 136 45 148 140 39 124 118 80 150 44 30 31 52 120 4])';
D.fn = @(xh) aluffi_pentini_orig(xh);
PB.f_star      = -0.35238607380003645;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 1;
PB.x0_idx      = [52 20];   % protocol starting point: Halton point 1002 (index representation)
PB.x0_alt_idx  = [88 24;14 33;135 23;61 26;97 30;24 19;117 29;43 31;79 20];   % 9 alternative starting points (Halton indices 1002 + 63 k)
PB.core_pidx   = 2;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('aluffi_pentini_k_2D: n = 2, pattern K, m = 151, f* = -0.3523860738 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('aluffi_pentini_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('aluffi_pentini_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('aluffi_pentini_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('aluffi_pentini_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('aluffi_pentini_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('aluffi_pentini_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = aluffi_pentini_orig(x)
%
% Purpose:
%
%    Function aluffi_pentini is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -0.3523
%    Global minima = [-1.0465 0]'
%    Local minima (two)
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
f = 0.25*x(1)^4-0.5*x(1)^2+0.1*x(1)+0.5*x(2)^2;
end
%
% End of aluffi_pentini.

