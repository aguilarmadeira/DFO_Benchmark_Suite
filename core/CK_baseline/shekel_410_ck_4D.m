function PB = shekel_410_ck_4D()
%SHEKEL_410_CK_4D  DFO Benchmark Suite v2 core: shekel_410, dimension 4, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -10.53640981669204 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K
%
lb_orig = [0 0 0 0];
ub_orig = [10 10 10 10];
lb_C_work = [0 0 0 0];
ub_C_work = [10 10 10 10];
T_break = [NaN 4 NaN 4];   % breakpoints of T (NaN for continuous)
T_image = [NaN 4.0005929374560214 NaN 3.9995098024822462];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'shekel_410_ck_4D';
PB.base_function = 'shekel_410';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of shekel_410 (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1];
PB.cardinality = [Inf 151 Inf 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shekel_410_ck_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([147 23 126 149 3 29 119 62 123 49 116 58 26 146 89 54 102 78 82 121 41 105 11 83 109 125 77 36 122 70 127 98 32 52 151 48 80 115 90 88 59 150 14 22 30 44 6 9 128 92 45 129 37 42 12 60 34 53 66 64 81 135 106 25 79 5 18 76 148 16 24 65 140 141 10 40 113 96 51 19 112 84 15 110 124 117 21 38 131 95 50 73 27 103 33 17 143 93 56 132 130 108 7 87 134 142 67 71 20 1 68 57 137 111 31 86 72 39 43 75 118 63 85 138 114 100 99 107 101 97 133 91 120 35 94 139 136 46 61 28 47 13 4 69 8 144 104 74 2 55 145])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([72 142 42 144 130 118 151 64 67 38 70 126 91 59 99 97 94 50 109 104 41 53 36 68 149 133 34 81 147 116 26 61 127 146 139 48 49 113 93 90 52 114 137 23 73 148 77 128 108 30 115 136 39 44 75 122 37 56 84 112 60 10 125 131 54 8 22 43 74 19 28 69 32 85 14 88 66 21 55 24 119 129 20 7 145 78 27 45 143 120 57 102 35 123 40 106 12 4 65 100 5 80 17 98 46 25 79 134 33 15 82 124 76 13 47 103 89 58 62 95 18 87 110 71 16 2 1 11 141 3 29 150 31 83 9 138 63 117 135 111 6 121 140 92 51 101 107 105 86 132 96])';
D.fn = @(xh) shekel_410_orig(xh);
PB.f_star      = -10.536409816692045;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [5.9423828125 133 9.7312000000000012 136];   % protocol starting point: Halton point 1049 (index representation)
PB.x0_alt_idx  = [1.0205078125 122 5.1712000000000007 142;9.1064453125 127 0.22720000000000001 149;4.1845703125 129 7.1872000000000007 135;6.5869140625 119 2.1664000000000003 141;1.6650390625 128 9.1264000000000003 148;7.8759765625 132 4.1824000000000003 130;2.9541015625 120 1.5423999999999998 137;5.3271484375 124 6.1984000000000004 143;0.4052734375 133 3.5583999999999993 150];   % 9 alternative starting points (Halton indices 1049 + 63 k)
PB.core_pidx   = 49;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shekel_410_ck_4D: n = 4, pattern CK, m = 151, f* = -10.53640982 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('shekel_410_ck_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('shekel_410_ck_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('shekel_410_ck_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
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
if ~iscell(xr) || numel(xr) ~= 4, error('shekel_410_ck_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_410_ck_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('shekel_410_ck_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shekel_410_ck_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shekel_410_ck_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = shekel_410_orig(x)
%
% Purpose:
%
%    Function shekel_4 is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value  
%    Global minima = 4*ones(4,1)
%    Local minima (m local minima at x = A(:,1:4))
%    Search domain: 0 <= x <= 10 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
%    Cases considered: m = 5 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 7 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 10 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
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
m = 10;
A = [4 4 4 4 0.1
    1 1 1 1 0.2
    8 8 8 8 0.2
    6 6 6 6 0.4
    3 7 3 7 0.4
    2 9 2 9 0.6
    5 5 3 3 0.3
    8 1 8 1 0.7
    6 2 6 2 0.5
    7 3.6 7 3.6 0.5];
f = -sum(1./(sum((repmat(x',m,1)-A(1:m,1:4)).^2,2)+A(1:m,5)),1);
end
%
% End of shekel_4.

