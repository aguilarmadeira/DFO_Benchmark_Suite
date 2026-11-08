function PB = kowalik_ck_4D()
%KOWALIK_CK_4D  DFO Benchmark Suite v2.1.0 core: kowalik, dimension 4, pattern CK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0.0003074859878056062 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K
%
lb_orig = [0.059979711914062504 0.081743057460752944 0 0.056734110787172012];
ub_orig = [0.41999999999999998 0.41999999999999998 0.34746458879999997 0.41999999999999998];
lb_C_work = [0.059979711914062483 0.081743057460752944 0 0.056734110787172012];
ub_C_work = [0.41999999999999998 0.41999999999999998 34.746458879999999 0.41999999999999998];
T_break = [NaN 0.189985279073312 NaN 0.13665260641399418];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0.1908362411132952 NaN 0.1357659904379078];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'kowalik_ck_4D';
PB.base_function = 'kowalik';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CK adaptation of kowalik (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100];
PB.cardinality = [Inf 151 Inf 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'kowalik_ck_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [0.19283345290075488 NaN 12.311730009947439 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [0.1928334529007549 NaN 0.1231173000994744 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([14 97 127 96 88 54 51 125 18 5 91 126 21 114 52 131 86 4 70 150 28 105 149 46 23 61 6 72 13 98 121 24 48 77 136 44 64 20 56 117 69 49 40 3 120 137 111 109 74 142 8 65 146 47 12 95 116 62 106 103 35 118 92 68 76 31 144 60 66 27 134 84 101 90 138 16 85 122 7 45 75 15 25 102 145 36 63 78 19 37 32 50 141 53 139 130 99 112 108 38 135 81 10 73 71 94 129 11 57 148 26 82 39 151 34 29 110 147 58 9 113 59 128 41 17 79 100 123 30 124 1 132 143 115 140 2 89 133 104 22 87 119 67 42 33 43 80 83 107 93 55])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([143 6 75 120 15 62 124 92 137 52 138 8 129 22 87 73 36 5 133 9 95 94 74 146 136 148 151 11 55 16 43 89 30 34 88 118 79 67 108 27 33 54 86 4 51 21 28 114 102 58 117 18 45 112 14 24 32 100 56 121 142 48 113 68 110 35 116 78 1 122 119 71 128 131 135 140 37 63 12 115 80 99 101 3 38 107 109 104 123 125 39 57 127 60 81 49 19 90 103 47 130 91 20 84 82 126 42 23 69 97 139 96 53 147 50 141 10 98 77 134 149 83 25 66 44 106 2 26 64 31 132 41 150 29 72 59 17 144 40 105 46 85 61 65 76 145 70 93 13 111 7])';
D.fn = @(xh) kowalik_orig(xh);
PB.f_star      = 0.00030748598780560617;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.10515803908109664 105 21.303748868505604 151];   % protocol starting point: Halton point 1028 (index representation)
PB.x0_alt_idx  = [0.33298337763547897 109 12.130683724185602 133;0.15156690434217451 114 28.308634978713599 140;0.24438463486433026 101 19.1355698343936 146;0.060858667695522291 111 1.6900677599232001 132;0.41525363878011701 115 25.873603140403201 139;0.23805615323781965 105 8.6949538701312008 142;0.32665489600896835 107 32.878489250611203 149;0.14945741046667099 116 15.699839980339199 134;0.36989952045679092 105 6.5267748360192002 141];   % 9 alternative starting points (Halton indices 1028 + 63 k)
PB.core_pidx   = 28;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('kowalik_ck_4D: n = 4, pattern CK, m = 151, f* = 0.0003074859878 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('kowalik_ck_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('kowalik_ck_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('kowalik_ck_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
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
if ~iscell(xr) || numel(xr) ~= 4, error('kowalik_ck_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('kowalik_ck_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('kowalik_ck_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('kowalik_ck_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('kowalik_ck_4D: x%d unknown label %s.', i, lab); end
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

