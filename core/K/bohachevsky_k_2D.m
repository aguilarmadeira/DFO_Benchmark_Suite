function PB = bohachevsky_k_2D()
%BOHACHEVSKY_K_2D  DFO Benchmark Suite v2 core: bohachevsky, dimension 2, pattern K, strategy baseline.
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
lb_orig = [-50 -50];
ub_orig = [50 50];
lb_C_work = [-50 -50];
ub_C_work = [50 50];
T_break = [0 0];   % breakpoints of T (NaN for continuous)
T_image = [0 0];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'bohachevsky_k_2D';
PB.base_function = 'bohachevsky';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of bohachevsky (F = f o T).';
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
D = struct('name', 'bohachevsky_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([105 55 25 50 93 64 140 143 23 54 76 127 11 148 72 82 65 150 151 149 68 99 10 75 36 77 79 120 119 43 102 95 52 12 108 16 21 117 109 33 69 146 97 5 9 118 94 103 29 61 71 122 101 67 73 42 83 134 132 17 35 41 74 89 84 13 110 116 37 124 136 114 34 115 123 145 111 113 57 4 112 90 47 2 142 3 104 32 125 39 129 20 130 87 7 98 1 51 40 60 18 131 81 22 62 14 44 106 30 92 19 24 38 15 45 88 28 100 49 141 63 139 66 58 53 59 85 86 8 46 126 6 78 135 70 137 121 80 96 133 26 27 107 31 91 56 128 48 144 147 138])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([139 12 57 149 126 51 104 84 143 106 8 117 142 92 135 82 75 3 124 112 32 9 134 78 83 115 63 26 145 89 44 123 127 148 66 137 125 29 49 105 38 52 58 80 141 70 101 87 85 61 46 1 19 21 88 62 48 2 7 18 121 43 74 54 122 131 77 140 39 150 65 17 118 47 69 10 144 40 119 103 132 90 93 138 30 6 111 35 56 42 95 24 107 91 11 110 5 55 113 128 146 130 86 28 68 20 50 4 36 98 147 31 45 22 53 97 120 109 59 37 73 108 79 71 67 72 99 102 25 64 15 27 100 94 96 13 16 116 133 34 76 81 33 114 151 23 129 136 14 60 41])';
D.fn = @(xh) bohachevsky_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [34 120];   % protocol starting point: Halton point 1004 (index representation)
PB.x0_alt_idx  = [125 124;51 133;88 123;15 126;135 130;61 119;98 129;24 131;116 120];   % 9 alternative starting points (Halton indices 1004 + 63 k)
PB.core_pidx   = 4;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('bohachevsky_k_2D: n = 2, pattern K, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('bohachevsky_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('bohachevsky_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('bohachevsky_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('bohachevsky_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('bohachevsky_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('bohachevsky_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = bohachevsky_orig(x)
%
% Purpose:
%
%    Function bohachevsky is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -50 <= x <= 50 (Montaz et al. (2005))
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
f = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2))+0.3;
end
%
% End of bohachevsky.

