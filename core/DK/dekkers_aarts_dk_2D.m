function PB = dekkers_aarts_dk_2D()
%DEKKERS_AARTS_DK_2D  DFO Benchmark Suite v2.1.0 core: dekkers_aarts, dimension 2, pattern DK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -24776.51834231769 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:K
%
lb_orig = [-20 -12.145811106030584];
ub_orig = [14.025146484375 20];
lb_C_work = [-20 -12.145811106030584];
ub_C_work = [14.025146484375 20];
T_break = [-0.038580729166667993 14.856670223035108];   % breakpoints of T (NaN for continuous)
T_image = [0 14.94511215764898];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'dekkers_aarts_dk_2D';
PB.base_function = 'dekkers_aarts';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: DK adaptation of dekkers_aarts (F = f o T).';
PB.nvar        = 2;
PB.nobj        = 1;
PB.var_types   = {'D','K'};
PB.type_counts = struct('C', 0, 'D', 1, 'K', 1);
PB.proportion  = [0 1 1];
PB.pattern     = 'DK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [151 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 151;
PB.domain_size = 22801;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'dekkers_aarts_dk_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([94 104 82 92 65 55 64 115 107 21 46 67 74 130 29 105 96 137 95 2 40 116 146 34 84 15 37 71 111 39 61 3 98 85 68 147 19 48 54 60 36 120 57 5 10 62 88 81 32 38 138 122 25 28 128 52 124 106 18 145 51 93 101 113 149 35 44 144 143 102 22 119 69 7 97 50 63 108 117 56 79 151 123 66 136 135 118 24 47 103 16 30 133 131 129 76 31 150 86 148 73 17 109 91 126 110 45 6 41 12 132 72 127 87 20 80 26 13 142 90 125 4 121 1 139 141 100 43 27 112 11 53 59 33 78 9 70 49 23 134 42 89 77 58 75 140 114 8 83 99 14])';
D.fn = @(xh) dekkers_aarts_orig(xh);
PB.f_star      = -24776.51834231769;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 1;
PB.x0_idx      = [48 109];   % protocol starting point: Halton point 1010 (index representation)
PB.x0_alt_idx  = [83 113;9 103;145 112;71 114;107 104;33 108;126 118;52 105;88 109];   % 9 alternative starting points (Halton indices 1010 + 63 k)
PB.core_pidx   = 10;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('dekkers_aarts_dk_2D: n = 2, pattern DK, m = 151, f* = -24776.51834 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('dekkers_aarts_dk_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.lb(1) + (k - 1) * (D.ub(1) - D.lb(1)) / (D.m - 1);
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
if ~iscell(xr) || numel(xr) ~= 2, error('dekkers_aarts_dk_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('dekkers_aarts_dk_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('dekkers_aarts_dk_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('dekkers_aarts_dk_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('dekkers_aarts_dk_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = dekkers_aarts_orig(x)
%
% Purpose:
%
%    Function dekkers_aarts is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -24771
%    Global minima = [0 15]', [0 -15]'
%    Local minima (three)
%    Search domain: -20 <= x <= 20 (Montaz et al. (2005))
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
f = 10^5*x(1)^2+x(2)^2-(x(1)^2+x(2)^2)^2+1/10^5*(x(1)^2+x(2)^2)^4;
end
%
% End of dekkers_aarts.

