function PB = kowalik_cdk_4D()
%KOWALIK_CDK_4D  DFO Benchmark Suite v2 core: kowalik, dimension 4, pattern CDK, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0.0003074859878056062 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C
%
lb_orig = [0 0 0 0];
ub_orig = [0.41999999999999998 0.41999999999999998 0.41999999999999998 0.41999999999999998];
lb_C_work = [0 0 0 0];
ub_C_work = [0.41999999999999998 0.41999999999999998 0.41999999999999998 42000000];
T_break = [NaN 0.19039999999999999 0.1232 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0.1908362411132952 0.1231173000994744 NaN];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
PB.name        = 'kowalik_cdk_4D';
PB.base_function = 'kowalik';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CDK adaptation of kowalik (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C'};
PB.type_counts = struct('C', 2, 'D', 1, 'K', 1);
PB.proportion  = [2 1 1];
PB.pattern     = 'CDK';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 100000000];
PB.cardinality = [Inf 151 151 Inf];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 151;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'kowalik_cdk_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([14 97 127 96 88 54 51 125 18 5 91 126 21 114 52 131 86 4 70 150 28 105 149 46 23 61 6 72 13 98 121 24 48 77 136 44 64 20 56 117 69 49 40 3 120 137 111 109 74 142 8 65 146 47 12 95 116 62 106 103 35 118 92 68 76 31 144 60 66 27 134 84 101 90 138 16 85 122 7 45 75 15 25 102 145 36 63 78 19 37 32 50 141 53 139 130 99 112 108 38 135 81 10 73 71 94 129 11 57 148 26 82 39 151 34 29 110 147 58 9 113 59 128 41 17 79 100 123 30 124 1 132 143 115 140 2 89 133 104 22 87 119 67 42 33 43 80 83 107 93 55])';
D.fn = @(xh) kowalik_orig(xh);
PB.f_star      = 0.00030748598780560617;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.052705078124999999 105 93 41912536.443148687];   % protocol starting point: Halton point 1028 (index representation)
PB.x0_alt_idx  = [0.318486328125 109 53 37032069.970845483;0.106845703125 114 123 38868804.664723031;0.21512695312499999 101 84 40705539.358600579;0.001025390625 111 8 36664723.032069974;0.41446289062499997 115 113 38501457.725947522;0.20774414062499999 105 39 39498542.274052471;0.31110351562499999 107 143 41335276.967930026;0.104384765625 116 69 37294460.641399413;0.36155273437499996 105 29 39131195.335276961];   % 9 alternative starting points (Halton indices 1028 + 63 k)
PB.core_pidx   = 28;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('kowalik_cdk_4D: n = 4, pattern CDK, m = 151, f* = 0.0003074859878 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('kowalik_cdk_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('kowalik_cdk_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('kowalik_cdk_4D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
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
if ~iscell(xr) || numel(xr) ~= 4, error('kowalik_cdk_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('kowalik_cdk_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('kowalik_cdk_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('kowalik_cdk_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('kowalik_cdk_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = kowalik_orig(x)
%
% Purpose:
%
%    Function kowalik is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 3.0748*10^-4
%    Global minima = [0.192 0.190 0.123 0.135]'
%    Local minima 
%    Search domain: 0 <= x <= 0.42 (Montaz et al. (2005))
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
a = [0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b = [0.25 0.5 1 2 4 6 8 10 12 14 16];
f = sum((a-x(1).*(1+x(2).*b)./(1+x(3).*b+x(4).*b.^2)).^2,2);
end
%
% End of kowalik.

