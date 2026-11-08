function PB = tenn_local_minima_k_2D()
%TENN_LOCAL_MINIMA_K_2D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 2, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 2.355817157679958e-31 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K
%
lb_orig = [-10 -10];
ub_orig = [6.6461669921875002 6.877569476197734];
lb_C_work = [-10 -10];
ub_C_work = [6.6461669921875002 6.877569476197734];
T_break = [0.98647021484375053 1.0266787244491873];   % breakpoints of T (NaN for continuous)
T_image = [1 1];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'tenn_local_minima_k_2D';
PB.base_function = 'tenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of tenn_local_minima (F = f o T).';
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
D = struct('name', 'tenn_local_minima_k_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([95 58 68 37 146 65 33 106 114 25 100 138 81 111 34 42 74 80 125 77 6 66 71 49 50 132 90 46 82 20 89 69 140 10 51 115 27 97 4 75 15 91 62 61 143 85 86 131 36 92 113 110 9 13 39 107 7 119 137 134 122 48 55 93 31 53 1 116 83 102 139 136 78 41 96 40 63 148 150 29 70 17 151 105 72 22 59 79 147 88 23 47 141 127 18 8 2 98 73 11 24 3 126 28 101 103 54 142 149 26 118 104 14 76 16 94 123 99 112 128 121 35 135 64 133 45 84 130 129 56 52 60 43 144 87 120 38 19 12 124 145 67 44 108 5 117 30 109 21 57 32])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([18 67 69 94 114 28 46 7 87 146 63 19 39 147 53 68 45 76 80 137 124 126 117 131 122 38 116 128 36 135 30 90 64 142 92 129 130 79 145 57 16 115 103 35 44 120 134 52 61 144 123 109 85 14 41 141 93 101 40 2 37 139 125 97 89 138 106 86 83 5 136 58 127 43 10 42 133 23 107 32 118 20 149 91 75 25 62 82 104 143 151 148 9 88 22 12 6 4 78 15 150 8 50 110 1 119 60 73 27 33 105 111 21 84 24 102 96 108 3 66 11 47 100 77 56 59 99 71 17 74 72 81 65 132 113 13 70 54 49 26 55 121 112 48 98 95 31 51 34 29 140])';
D.fn = @(xh) tenn_local_minima_orig(xh);
PB.f_star      = 2.3558171576799584e-31;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [78 69];   % protocol starting point: Halton point 1057 (index representation)
PB.x0_alt_idx  = [5 78;147 82;73 70;109 74;35 84;129 73;55 75;90 80;16 69];   % 9 alternative starting points (Halton indices 1057 + 63 k)
PB.core_pidx   = 57;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('tenn_local_minima_k_2D: n = 2, pattern K, m = 151, f* = 2.355817158e-31 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('tenn_local_minima_k_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('tenn_local_minima_k_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_k_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_k_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('tenn_local_minima_k_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('tenn_local_minima_k_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = tenn_local_minima_orig(x)
%
% Purpose:
%
%    Function tenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (10^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
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
dim = size(x,1);
f   = pi/dim*sum(((x(1:dim-1)-1).^2).*(1+10*sin(pi.*x(2:dim)).^2),1);
f   = f + pi/dim*(10*sin(pi*x(1))^2+(x(dim)-1)^2);
end
%
% End of tenn_local_minima.

