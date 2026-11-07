function PB = miele_cantrell_ck_4D()
%MIELE_CANTRELL_CK_4D  DFO Benchmark Suite v2 core: miele_cantrell, dimension 4, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K
%
lb_orig = [-10 -10 -10 -10];
ub_orig = [10 10 10 10];
lb_C_work = [-10 -10 -10 -10];
ub_C_work = [10 10 10 10];
T_break = [NaN 1.0666666666666667 NaN 1.0666666666666667];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 NaN 1];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'miele_cantrell_ck_4D';
PB.base_function = 'miele_cantrell';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of miele_cantrell (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1];
PB.cardinality = [Inf 151 Inf 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'miele_cantrell_ck_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([77 122 104 101 150 18 6 128 114 56 57 68 102 78 82 129 36 48 71 19 107 52 111 8 79 123 141 26 127 3 96 110 15 73 66 131 86 121 140 92 16 49 97 31 119 4 11 14 70 139 44 81 134 142 105 9 76 94 41 58 25 137 133 80 7 23 74 53 126 22 28 63 116 148 90 17 34 62 115 51 145 24 35 10 47 124 84 108 42 67 65 85 32 93 20 87 132 37 72 1 106 50 83 138 54 89 27 12 60 99 45 147 151 98 88 95 64 55 38 5 120 112 69 118 130 59 135 39 117 113 43 61 40 29 100 143 125 136 149 91 21 33 46 103 2 13 109 144 30 146 75])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([30 85 131 64 146 41 51 62 84 59 50 37 144 129 118 139 70 122 36 142 55 87 5 8 47 17 22 27 110 137 11 145 48 76 42 97 140 106 148 68 88 89 134 33 133 6 13 79 149 75 98 4 150 147 132 116 130 35 44 1 96 72 43 54 100 151 38 80 63 120 32 67 108 77 128 21 136 66 16 56 123 29 143 125 53 61 3 18 49 91 109 92 107 7 28 95 138 46 117 10 99 60 93 119 65 101 39 24 73 127 58 12 90 113 104 112 81 74 57 25 124 14 94 20 135 86 2 71 23 19 82 102 83 78 15 103 52 69 105 34 121 141 40 126 45 115 111 9 114 31 26])';
D.fn = @(xh) miele_cantrell_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-8.740234375 38 -0.93759999999999977 65];   % protocol starting point: Halton point 1032 (index representation)
PB.x0_alt_idx  = [7.666015625 43 -6.2176000000000009 72;-2.412109375 47 3.094400000000002 79;2.744140625 35 -2.1856 85;-7.451171875 44 6.9727999999999994 71;5.205078125 48 1.6928000000000019 77;-4.873046875 38 -8.1951999999999998 81;0.283203125 40 5.7248000000000037 66;-9.970703125 49 -4.1631999999999998 73;9.716796875 39 9.7568000000000019 80];   % 9 alternative starting points (Halton indices 1032 + 63 k)
PB.core_pidx   = 32;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('miele_cantrell_ck_4D: n = 4, pattern CK, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('miele_cantrell_ck_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('miele_cantrell_ck_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('miele_cantrell_ck_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
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
if ~iscell(xr) || numel(xr) ~= 4, error('miele_cantrell_ck_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('miele_cantrell_ck_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('miele_cantrell_ck_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('miele_cantrell_ck_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('miele_cantrell_ck_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = miele_cantrell_orig(x)
%
% Purpose:
%
%    Function miele_cantrell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [0 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = (exp(x(1))-x(2))^4 + 100*(x(2)-x(3))^6 + tan(x(3)-x(4))^4 + x(1)^8;
end
%
% End of miele_cantrell.

