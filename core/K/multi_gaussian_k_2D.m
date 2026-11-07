function PB = multi_gaussian_k_2D()
%MULTI_GAUSSIAN_K_2D  DFO Benchmark Suite v2.1.0 core: multi_gaussian, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -1.296954045953779 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-1.3062646484375 -1.3355687649240462];
ub_orig = [2 2];
lb_C_work = [-1.3062646484375 -1.3355687649240462];
ub_C_work = [2 2];
T_break = [-0.0058005533854168512 -0.023578384053921431];   % breakpoints of T (NaN for continuous)
T_image = [-0.01354066395522768 -0.01354066351404425];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'multi_gaussian_k_2D';
PB.base_function = 'multi_gaussian';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of multi_gaussian (F = f o T).';
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
D = struct('name', 'multi_gaussian_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([101 72 135 132 39 110 17 55 146 54 12 33 49 23 64 73 58 126 47 104 90 5 105 82 32 120 149 100 127 7 117 15 119 151 57 79 59 62 145 16 121 28 134 50 116 11 18 87 141 112 1 118 38 42 99 65 27 111 107 75 68 46 77 63 103 80 48 136 4 138 123 6 128 95 9 34 52 91 30 150 96 114 56 31 69 133 125 115 67 3 2 21 61 106 51 22 93 71 8 37 25 88 19 86 97 142 74 60 109 41 98 53 20 139 144 45 14 143 113 129 81 83 43 94 147 40 85 35 122 124 66 89 78 84 26 131 102 140 44 137 130 29 10 148 13 108 92 36 76 70 24])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([27 69 131 65 78 142 40 127 122 71 18 143 62 76 129 46 31 77 96 125 119 42 85 61 15 82 84 148 5 7 14 16 56 104 57 67 128 114 130 17 149 29 32 50 1 133 19 22 101 59 2 112 88 150 9 118 113 115 63 98 108 110 74 102 21 89 49 68 90 38 45 8 37 53 11 36 116 4 94 73 146 47 58 138 72 139 111 141 70 121 132 25 64 33 123 26 97 75 12 41 30 144 24 92 103 140 80 66 6 48 105 60 28 135 44 126 136 20 3 151 91 95 55 145 106 54 100 51 10 13 86 137 99 109 52 23 134 34 81 35 39 83 79 120 107 87 93 124 117 43 147])';
D.fn = @(xh) multi_gaussian_orig(xh);
PB.f_star      = -1.2969540459537789;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [85 88];   % protocol starting point: Halton point 1033 (index representation)
PB.x0_alt_idx  = [12 93;133 97;59 85;95 94;21 98;114 88;41 90;76 99;2 89];   % 9 alternative starting points (Halton indices 1033 + 63 k)
PB.core_pidx   = 33;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('multi_gaussian_k_2D: n = 2, pattern K, m = 151, f* = -1.296954046 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('multi_gaussian_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('multi_gaussian_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('multi_gaussian_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('multi_gaussian_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('multi_gaussian_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('multi_gaussian_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = multi_gaussian_orig(x)
%
% Purpose:
%
%    Function multi_gaussian is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -1.29695
%    Global minima = [-0.01356 -0.01356]'
%    Local minima (five)
%    Search domain: -2 <= x <= 2 (Montaz et al. (2005))
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
A =[0.5 0 0 0.1
    1.2 1 0 0.5
    1 0 -0.5 0.5
    1 -0.5 0 0.5
    1.2 0 1 0.5];
f = -sum(A(:,1).*exp(-((x(1)-A(:,2)).^2+(x(2)-A(:,3)).^2)./A(:,4).^2),1);
end
%
% End of multi_gaussian.

