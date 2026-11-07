function PB = branin_hoo_dk_2D()
%BRANIN_HOO_DK_2D  DFO Benchmark Suite v2.1.0 core: branin_hoo, dimension 2, pattern DK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0.3978873577297382 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:K
%
lb_orig = [-5 2.1286541685718641];
ub_orig = [7.2836486816406252 15];
lb_C_work = [-5 2.1286541685718641];
ub_C_work = [7.2836486816406252 15];
T_break = [-3.1165072021484379 12.254112889295332];   % breakpoints of T (NaN for continuous)
T_image = [-3.1415926535897931 12.275];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'branin_hoo_dk_2D';
PB.base_function = 'branin_hoo';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: DK adaptation of branin_hoo (F = f o T).';
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
D = struct('name', 'branin_hoo_dk_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([110 79 44 135 119 53 17 88 144 8 145 21 36 127 80 77 128 147 96 43 87 32 14 16 38 23 61 2 51 70 71 132 15 50 20 22 107 35 136 40 57 7 151 108 140 26 62 117 134 141 98 10 66 148 46 69 3 138 19 30 68 92 78 4 18 76 37 81 83 6 105 34 54 101 12 89 41 60 1 28 49 64 104 29 116 33 52 63 84 106 146 56 82 31 121 109 42 93 85 102 123 120 133 73 11 72 103 131 122 129 86 142 111 90 9 55 113 74 24 5 139 91 143 125 59 149 99 114 47 94 58 75 27 39 48 126 25 13 45 118 95 130 65 124 137 97 100 112 67 150 115])';
D.fn = @(xh) branin_hoo_orig(xh);
PB.f_star      = 0.39788735772973816;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 2;
PB.x0_idx      = [109 36];   % protocol starting point: Halton point 1005 (index representation)
PB.x0_alt_idx  = [32 41;126 50;51 40;90 42;13 47;136 36;61 45;99 47;22 37];   % 9 alternative starting points (Halton indices 1005 + 63 k)
PB.core_pidx   = 5;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('branin_hoo_dk_2D: n = 2, pattern DK, m = 151, f* = 0.3978873577 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('branin_hoo_dk_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('branin_hoo_dk_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('branin_hoo_dk_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('branin_hoo_dk_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('branin_hoo_dk_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('branin_hoo_dk_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = branin_hoo_orig(x)
%
% Purpose:
%
%    Function brannin_hoo is the function described in
%    Garcia-Palomares(2012).
%
%    dim = 2
%    Minimum global value =  5/(4*pi)
%    Global minima = ([-pi 12.275]', [pi 2.275]', [3*pi 2.475]')
%    Local minima (equal to global minima)
%    Search domain: -5 <= x <= 20 (Garcia-Palomares(2012))
%
%                    -5 <= x(1) <= 10 
%                     0 <= x(2) <= 15 (Huyer and Neumaier (1999))
%                                     (Huyer and Neumaier (2008))
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
f = (x(2)-1.275*(x(1)/pi)^2+5*x(1)/pi-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end
%
% End of branin_hoo.

