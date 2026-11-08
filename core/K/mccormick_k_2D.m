function PB = mccormick_k_2D()
%MCCORMICK_K_2D  DFO Benchmark Suite v2.1.0 core: mccormick, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -1.913222954981037 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-1.5 -3];
ub_orig = [3.1183013916015625 2.0500198140527357];
lb_C_work = [-1.5 -3];
ub_C_work = [3.1183013916015625 2.0500198140527357];
T_break = [-0.54555104573567703 -1.5523276533048824];   % breakpoints of T (NaN for continuous)
T_image = [-0.54719755371679324 -1.547197552675593];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'mccormick_k_2D';
PB.base_function = 'mccormick';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of mccormick (F = f o T).';
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
D = struct('name', 'mccormick_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([26 78 130 131 2 79 100 48 99 137 144 123 83 138 113 51 62 121 116 14 44 92 111 133 39 103 81 17 40 24 38 128 108 9 105 122 147 73 60 37 67 134 10 102 107 132 41 45 90 143 25 34 77 11 31 141 50 8 16 33 66 7 96 94 88 59 3 118 112 82 65 87 119 71 28 109 70 19 47 6 36 58 125 95 148 110 27 42 72 4 91 18 149 142 43 12 80 104 106 145 101 150 93 126 120 86 23 124 52 127 29 89 49 69 53 54 20 117 98 129 46 13 84 5 114 146 151 1 57 135 68 76 32 115 35 64 75 139 30 15 74 22 63 136 140 55 21 97 56 61 85])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([48 37 89 121 20 41 15 122 64 145 126 10 27 32 4 6 80 51 43 36 112 127 77 148 55 81 50 141 78 88 1 86 79 95 12 14 68 31 102 151 99 53 119 70 24 129 42 97 71 114 94 111 109 105 150 65 52 116 17 146 115 9 13 23 44 85 101 147 117 59 120 137 113 73 143 130 83 106 49 118 39 93 72 46 82 2 30 45 75 7 125 21 103 26 47 128 84 104 110 19 134 138 98 60 63 136 29 11 58 91 35 123 57 142 135 62 28 8 108 139 140 25 132 18 3 149 100 16 74 133 87 96 54 5 61 92 107 66 67 56 124 144 131 76 90 34 22 40 69 33 38])';
D.fn = @(xh) mccormick_orig(xh);
PB.f_star      = -1.9132229549810367;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [57 72];   % protocol starting point: Halton point 1030 (index representation)
PB.x0_alt_idx  = [96 76;20 80;115 68;39 77;78 82;2 71;150 73;76 83;111 72];   % 9 alternative starting points (Halton indices 1030 + 63 k)
PB.core_pidx   = 30;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('mccormick_k_2D: n = 2, pattern K, m = 151, f* = -1.913222955 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('mccormick_k_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
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
if ~iscell(xr) || numel(xr) ~= 2, error('mccormick_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('mccormick_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('mccormick_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('mccormick_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('mccormick_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = mccormick_orig(x)
%
% Purpose:
%
%    Function mccormick is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -1.9133
%    Global minima = [-0.547 -1.547]'
%    Local minima (two)
%    Search domain: -1.5 <= x(1) <= 4
%                     -3 <= x(2) <= 3 (Montaz et al. (2005))
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
f = sin(x(1)+x(2))+(x(1)-x(2))^2-3/2*x(1)+5/2*x(2)+1;
end
%
% End of mccormick.

