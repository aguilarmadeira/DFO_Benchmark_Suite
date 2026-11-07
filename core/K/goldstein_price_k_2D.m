function PB = goldstein_price_k_2D()
%GOLDSTEIN_PRICE_K_2D  DFO Benchmark Suite v2 core: goldstein_price, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 3 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-2 -2];
ub_orig = [2 2];
lb_C_work = [-2 -2];
ub_C_work = [2 2];
T_break = [0 -0.98666666666666669];   % breakpoints of T (NaN for continuous)
T_image = [0 -1];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'goldstein_price_k_2D';
PB.base_function = 'goldstein_price';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of goldstein_price (F = f o T).';
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
D = struct('name', 'goldstein_price_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([13 35 137 106 30 32 67 44 140 18 48 71 151 22 41 3 90 60 108 102 129 136 123 40 113 78 66 134 12 38 58 23 133 93 34 117 142 15 125 95 49 76 126 39 74 10 72 63 53 33 143 2 91 130 47 94 92 145 147 131 149 14 68 114 107 121 4 19 139 148 88 81 105 56 9 146 87 70 17 150 116 16 101 96 124 21 122 59 86 5 83 98 120 1 20 77 111 8 25 50 138 80 6 64 43 45 112 57 99 104 42 118 24 28 103 79 36 84 54 97 61 46 82 109 55 100 85 26 119 29 73 37 27 75 62 135 132 128 127 69 52 7 65 141 89 31 110 144 115 51 11])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([139 64 37 126 70 145 4 149 141 39 72 132 44 43 41 52 94 76 121 118 92 14 17 55 112 60 67 10 98 15 74 123 77 113 35 110 129 148 106 125 138 32 50 40 56 75 18 108 53 34 65 3 104 59 109 2 151 22 101 78 8 16 20 83 130 116 6 21 136 84 142 79 117 124 103 33 135 73 82 127 96 19 134 144 28 24 115 62 89 146 86 12 97 5 150 81 46 13 29 54 30 85 11 69 48 51 90 119 105 1 143 91 31 36 140 88 45 131 63 107 71 58 95 122 68 114 100 42 7 47 93 57 49 99 87 137 27 26 25 147 102 61 133 80 38 9 111 66 128 23 120])';
D.fn = @(xh) goldstein_price_orig(xh);
PB.f_star      = 3;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [113 81];   % protocol starting point: Halton point 1021 (index representation)
PB.x0_alt_idx  = [36 70;131 75;56 84;94 72;18 76;141 80;65 68;104 77;27 81];   % 9 alternative starting points (Halton indices 1021 + 63 k)
PB.core_pidx   = 21;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('goldstein_price_k_2D: n = 2, pattern K, m = 151, f* = 3 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('goldstein_price_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('goldstein_price_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('goldstein_price_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('goldstein_price_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('goldstein_price_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('goldstein_price_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = goldstein_price_orig(x)
%
% Purpose:
%
%    Function goldstein_price is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = 3
%    Global minima = [0 -1]'
%    Local minima (four)
%    Search domain: -2 <= x <= 2 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
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
f = 1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2);
f = f * (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2));
end
%
% End of goldstein_price.

