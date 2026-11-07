function PB = threehumpcamel_k_2D()
%THREEHUMPCAMEL_K_2D  DFO Benchmark Suite v2 core: threehumpcamel, dimension 2, pattern K, strategy baseline.
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
lb_orig = [-5 -5];
ub_orig = [5 5];
lb_C_work = [-5 -5];
ub_C_work = [5 5];
T_break = [0 0];   % breakpoints of T (NaN for continuous)
T_image = [0 0];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'threehumpcamel_k_2D';
PB.base_function = 'threehumpcamel';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: K adaptation of threehumpcamel (F = f o T).';
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
D = struct('name', 'threehumpcamel_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([32 37 9 14 81 139 124 129 90 55 147 67 89 132 3 4 47 30 142 58 40 112 128 29 149 16 134 91 5 110 115 46 79 97 99 85 12 61 63 34 104 151 57 100 136 116 141 10 21 94 113 11 48 96 64 76 72 103 101 49 69 38 66 59 138 130 95 17 144 108 111 52 19 39 121 83 65 109 62 105 78 126 86 36 2 150 80 125 54 27 56 82 43 135 60 106 107 140 35 137 119 118 133 1 145 146 120 68 148 42 114 20 50 6 70 44 41 74 131 51 7 88 33 18 13 26 84 102 87 45 92 122 15 22 31 28 73 143 77 53 75 123 23 25 127 98 71 24 93 117 8])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([19 145 40 14 77 97 144 70 117 80 65 137 75 133 36 4 7 149 33 15 98 73 2 129 48 88 23 94 85 47 68 87 25 90 41 3 111 61 130 34 1 9 116 106 82 45 71 11 22 32 119 12 110 147 59 108 53 121 28 51 122 84 66 91 24 93 120 18 151 69 8 54 21 42 132 43 150 141 140 38 81 29 143 127 5 44 128 16 101 124 60 86 123 35 64 107 55 113 126 115 96 79 20 6 89 103 136 74 142 49 139 27 138 13 78 52 50 83 102 62 17 100 134 31 26 39 99 118 104 63 148 10 37 46 56 57 105 58 114 95 125 30 92 109 146 76 72 112 135 67 131])';
D.fn = @(xh) threehumpcamel_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [97 135];   % protocol starting point: Halton point 1061 (index representation)
PB.x0_alt_idx  = [23 144;117 149;43 137;79 141;5 150;147 140;73 142;109 146;35 136];   % 9 alternative starting points (Halton indices 1061 + 63 k)
PB.core_pidx   = 61;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('threehumpcamel_k_2D: n = 2, pattern K, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('threehumpcamel_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('threehumpcamel_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('threehumpcamel_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('threehumpcamel_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('threehumpcamel_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('threehumpcamel_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = threehumpcamel_orig(x)
%
% Purpose:
%
%    Function threehumpcamel is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima (three)
%    Search domain: -5 <= x <= 5 (Montaz et al. (2005))
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
f = 2*x(1)^2-1.05*x(1)^4+1/6*x(1)^6+x(1)*x(2)+x(2)^2;
end
%
% End of threehumpcamel.

