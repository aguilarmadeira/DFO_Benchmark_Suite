function PB = rosenbrock_k_2D()
%ROSENBROCK_K_2D  DFO Benchmark Suite v2.1.0 core: rosenbrock, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-5.1200000000000001 -5.1200000000000001];
ub_orig = [3.4252375000000002 3.6314412233907434];
lb_C_work = [-5.1200000000000001 -5.1200000000000001];
ub_C_work = [3.4252375000000002 3.6314412233907434];
T_break = [0.97560275000000019 1.0060088563735201];   % breakpoints of T (NaN for continuous)
T_image = [1 1];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'rosenbrock_k_2D';
PB.base_function = 'rosenbrock';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of rosenbrock (F = f o T).';
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
D = struct('name', 'rosenbrock_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([83 99 111 123 88 28 92 44 77 45 97 142 68 72 55 12 54 90 113 65 4 138 115 13 96 22 76 5 9 11 19 21 82 110 47 64 53 66 51 60 125 41 84 150 144 69 35 118 120 109 23 87 148 85 37 49 124 106 27 46 121 130 149 43 2 133 145 52 101 89 137 71 128 141 78 36 56 8 117 3 129 103 75 57 14 31 105 98 25 100 70 91 42 16 40 114 95 132 30 116 59 26 112 139 63 140 122 62 7 143 24 108 93 17 29 61 39 79 73 74 94 67 146 1 151 126 86 127 38 102 147 80 20 135 48 81 33 6 32 131 134 136 119 15 104 58 18 10 107 50 34])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([51 116 145 124 122 144 110 95 39 121 114 117 32 52 104 30 64 85 19 82 4 7 59 136 18 22 150 5 9 12 77 129 26 102 35 141 118 87 23 93 103 41 57 56 108 62 36 81 49 8 106 113 75 128 76 74 97 94 29 125 17 89 91 45 92 96 71 54 61 27 86 1 140 50 11 131 58 10 133 6 40 90 78 60 16 99 127 119 28 137 107 112 46 147 44 21 2 130 34 66 63 111 67 115 139 24 151 68 13 80 31 3 101 25 146 123 47 88 83 84 105 79 70 126 142 14 132 15 55 120 143 100 138 72 73 109 65 42 69 33 43 48 38 148 53 37 20 98 149 134 135])';
D.fn = @(xh) rosenbrock_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [81 44];   % protocol starting point: Halton point 1041 (index representation)
PB.x0_alt_idx  = [7 48;142 38;68 40;104 50;31 39;124 43;50 46;86 35;12 44];   % 9 alternative starting points (Halton indices 1041 + 63 k)
PB.core_pidx   = 41;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('rosenbrock_k_2D: n = 2, pattern K, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('rosenbrock_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('rosenbrock_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('rosenbrock_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('rosenbrock_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('rosenbrock_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('rosenbrock_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = rosenbrock_orig(x)
%
% Purpose:
%
%    Function rosenbrock is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (equal to global minima)
%    Search domain: -1000 <= x <= 1000 (Brachetti et al. (1997))
%                   -2.048 <= x <= 2.048 (Storn and Price (1997))
%                   -5.12 <= x <= 5.12 (Huyer and Neumaier (2008))
%    Cases considered: n = 2 Huyer and Neumaier (2008)
%                            Brachetti et al. (1997)
%                      n = 10 Garcia-Palomares  et al. (2006)
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
n    = size(x,1);
xaux = x(2:n);
f    = sum((x(1:n-1) - 1).^2 + 100.*(x(1:n-1).^2 - xaux).^2,1);
end
%
% End of rosenbrock.

