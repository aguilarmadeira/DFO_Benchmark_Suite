function PB = sixhumpcamel_dk_2D()
%SIXHUMPCAMEL_DK_2D  DFO Benchmark Suite v2.1.0 core: sixhumpcamel, dimension 2, pattern DK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -1.031628453488552 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:K
%
lb_orig = [-2.03322509765625 -1.3236620433876949];
ub_orig = [3 2];
lb_C_work = [-2.03322509765625 -1.3236620433876949];
ub_C_work = [3 2];
T_break = [0.08072944335937482 -0.70324512862199184];   % breakpoints of T (NaN for continuous)
T_image = [0.089842000000000005 -0.71265599999999996];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'sixhumpcamel_dk_2D';
PB.base_function = 'sixhumpcamel';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: DK adaptation of sixhumpcamel (F = f o T).';
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
D = struct('name', 'sixhumpcamel_dk_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([103 143 53 49 141 87 27 42 142 15 71 20 58 28 64 73 55 128 149 61 122 130 3 72 12 31 65 80 59 8 126 133 91 125 34 1 62 45 132 30 54 100 77 95 113 116 46 147 105 22 39 51 136 84 119 23 98 52 5 70 109 10 123 75 69 11 25 102 56 26 135 7 137 63 134 106 94 120 110 47 6 33 92 89 85 32 68 88 38 121 151 48 117 129 14 104 90 138 57 93 144 83 21 2 66 60 9 18 79 76 36 139 50 108 41 111 101 97 107 40 96 127 17 112 78 131 99 37 145 148 74 67 19 29 124 82 146 86 35 44 140 114 16 150 81 43 118 115 13 4 24])';
D.fn = @(xh) sixhumpcamel_orig(xh);
PB.f_star      = -1.0316284534885518;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 2;
PB.x0_idx      = [71 52];   % protocol starting point: Halton point 1054 (index representation)
PB.x0_alt_idx  = [110 61;34 65;129 53;53 58;92 67;16 56;139 59;62 63;101 52];   % 9 alternative starting points (Halton indices 1054 + 63 k)
PB.core_pidx   = 54;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('sixhumpcamel_dk_2D: n = 2, pattern DK, m = 151, f* = -1.031628453 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('sixhumpcamel_dk_2D: z must be a finite numeric vector with 2 entries.'); end
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
if ~iscell(xr) || numel(xr) ~= 2, error('sixhumpcamel_dk_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('sixhumpcamel_dk_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('sixhumpcamel_dk_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('sixhumpcamel_dk_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('sixhumpcamel_dk_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = sixhumpcamel_orig(x)
%
% Purpose:
%
%    Function sixhumpcamel is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = -1.0316285
%    Global minima: [-0.089842 0.712656]', [0.089842 -0.712656]'
%    Local minima (six)
%
%    Search domain: -2.5 <= x(1) <= 2.5 
%                   -1.5 <= x(2) <= 1.5 (Brachetti et al. (1997))
%
%    or
%
%                   -3 <= x(1) <= 3 
%                   -2 <= x(2) <= 2 (Huyer and Neumaier (1999))
%                                   (Huyer and Neumaier (2008))
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
f = (4 - 2.1*x(1)^2 + x(1)^4/3)*x(1)^2 + x(1)*x(2) + (4*x(2)^2 - 4)*x(2)^2;
end
%
% End of sixhumpcamel.

