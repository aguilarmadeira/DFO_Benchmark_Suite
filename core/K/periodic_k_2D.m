function PB = periodic_k_2D()
%PERIODIC_K_2D  DFO Benchmark Suite v2.1.0 core: periodic, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0.9 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-10 -10];
ub_orig = [6.3563232421875 6.0556216023980092];
lb_C_work = [-10 -10];
ub_C_work = [6.3563232421875 6.0556216023980092];
T_break = [0.031878255208333428 -0.045514606513233247];   % breakpoints of T (NaN for continuous)
T_image = [0 0];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'periodic_k_2D';
PB.base_function = 'periodic';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of periodic (F = f o T).';
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
D = struct('name', 'periodic_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([29 55 117 132 31 140 43 25 131 78 18 8 11 27 115 52 73 108 75 86 76 145 126 97 4 41 46 17 21 96 1 67 150 65 49 15 141 34 51 72 7 13 151 23 60 125 109 101 102 22 136 47 139 69 74 116 122 148 39 147 144 91 114 82 10 28 40 118 83 128 38 146 33 134 9 120 54 95 35 77 30 53 66 44 84 12 42 58 2 32 36 56 5 70 119 100 45 26 68 99 103 110 92 59 123 24 71 61 138 130 149 129 57 87 90 111 127 104 93 79 88 14 112 48 63 20 85 81 142 113 107 19 50 143 105 37 137 3 106 62 64 80 133 98 135 89 121 94 6 124 16])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([69 22 44 19 112 71 49 90 109 120 28 122 126 140 143 53 6 33 151 136 91 46 20 30 145 92 13 18 21 76 147 131 74 75 108 16 24 94 55 85 8 14 17 101 115 63 150 107 57 142 102 149 25 66 51 78 105 70 41 59 27 138 43 54 111 106 42 88 83 81 40 116 35 52 11 80 82 7 144 79 137 56 146 47 141 15 45 61 5 36 132 124 9 99 2 133 50 31 73 3 119 62 128 110 86 29 125 67 10 118 12 87 64 95 121 139 93 113 103 89 98 26 135 60 77 34 100 129 97 1 127 39 72 23 130 65 32 38 4 104 114 134 58 148 84 68 123 37 48 96 117])';
D.fn = @(xh) periodic_orig(xh);
PB.f_star      = 0.90000000000000002;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [104 111];   % protocol starting point: Halton point 1037 (index representation)
PB.x0_alt_idx  = [30 115;124 104;50 107;86 116;12 106;133 110;59 112;95 102;21 111];   % 9 alternative starting points (Halton indices 1037 + 63 k)
PB.core_pidx   = 37;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('periodic_k_2D: n = 2, pattern K, m = 151, f* = 0.9 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('periodic_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('periodic_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('periodic_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('periodic_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('periodic_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('periodic_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = periodic_orig(x)
%
% Purpose:
%
%    Function periodic is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0.9
%    Global minima = [0 0]'
%    Local minima (49 with value 1)
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
f = 1+sin(x(1))^2+sin(x(2))^2-0.1*exp(-x(1)^2-x(2)^2);
end
%
% End of periodic.

