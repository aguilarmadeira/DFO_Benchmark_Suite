function PB = shekel_foxholes_ck_5D()
%SHEKEL_FOXHOLES_CK_5D  DFO Benchmark Suite v2 core: shekel_foxholes, dimension 5, pattern CK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -10.40561723899244 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C
%
lb_orig = [0 0 0 0 0];
ub_orig = [10 10 10 10 10];
lb_C_work = [0 0 0 0 0];
ub_C_work = [10 10 1000 10 0.10000000000000001];
T_break = [NaN 9.1333333333333329 NaN 7.5999999999999996 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 9.1517257530637934 NaN 7.6208609576192838 NaN];   % their images
m = 151;
ProblemData = cell(1, 5);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
PB.name        = 'shekel_foxholes_ck_5D';
PB.base_function = 'shekel_foxholes';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of shekel_foxholes (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C'};
PB.type_counts = struct('C', 3, 'D', 0, 'K', 2);
PB.proportion  = [3 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01];
PB.cardinality = [Inf 151 Inf 151 Inf];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shekel_foxholes_ck_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([126 91 68 28 56 85 58 76 124 84 2 88 62 30 140 23 100 64 26 73 49 94 53 42 105 103 139 35 87 43 120 70 77 146 148 15 138 57 1 5 11 18 21 109 50 151 149 108 69 93 118 150 133 130 22 143 116 7 16 78 9 20 27 39 74 19 98 55 121 32 80 10 129 106 145 59 4 37 72 40 79 24 101 25 123 92 48 66 113 142 81 125 110 97 67 52 41 131 96 117 136 12 51 29 83 134 13 111 147 60 17 128 14 141 102 45 31 3 132 44 90 135 71 34 107 95 99 119 122 127 33 137 47 114 46 38 8 82 112 115 65 63 36 86 104 61 75 54 89 144 6])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([55 133 53 41 63 51 79 8 34 17 83 93 95 89 103 113 32 114 4 46 119 11 52 57 150 65 68 35 146 102 71 127 126 10 91 101 26 56 2 99 13 92 116 148 73 7 16 124 82 104 5 15 97 87 107 38 105 106 86 100 12 144 143 147 81 23 37 140 66 36 44 138 50 122 130 120 125 98 85 45 115 108 43 121 76 19 54 131 21 70 109 20 67 117 75 151 132 42 6 40 47 112 61 39 136 77 141 142 88 74 31 118 129 80 145 62 149 128 27 64 111 69 135 59 9 123 1 22 134 33 137 49 90 28 96 94 72 14 48 58 25 29 30 110 3 18 78 84 139 24 60])';
D.fn = @(xh) shekel_foxholes_orig(xh);
PB.f_star      = -10.405617238992441;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [3.4423828125 49 21.120000000000001 10 0.051840721262208879];   % protocol starting point: Halton point 1050 (index representation)
PB.x0_alt_idx  = [6.0205078125 39 717.12000000000012 17 0.020510894064613074;0.9814453125 43 222.72 2 0.097370398196844496;9.1845703125 46 918.72000000000003 9 0.066040570999248691;4.0869140625 35 416.63999999999999 16 0.043726521412471825;6.6650390625 45 152.63999999999999 22 0.011577078068437949;1.6259765625 49 618.24000000000012 5 0.088436582200669356;7.9541015625 36 354.23999999999995 11 0.057106755003073564;2.8271484375 41 819.84000000000003 18 0.034792705416296692;5.4052734375 50 555.84000000000003 3 0.0034628782187009089];   % 9 alternative starting points (Halton indices 1050 + 63 k)
PB.core_pidx   = 50;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shekel_foxholes_ck_5D: n = 5, pattern CK, m = 151, f* = -10.40561724 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('shekel_foxholes_ck_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('shekel_foxholes_ck_5D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('shekel_foxholes_ck_5D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('shekel_foxholes_ck_5D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
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
xr = cell(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 5, error('shekel_foxholes_ck_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_foxholes_ck_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('shekel_foxholes_ck_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shekel_foxholes_ck_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shekel_foxholes_ck_5D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = shekel_foxholes_orig(x)
%
% Purpose:
%
%    Function shekel_foxholes is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value
%    Global minima:  [8.025 9.152 5.114 7.621 4.564]', for n =5
%                    [8.025 9.152 5.114 7.621 4.564 4.771...
%                    2.996 6.126 0.734 4.982]' , for n = 10
%    Local minima (at least, twelve)
%    Search domain: 0 <= x <= 10 (Huyer and Neumaier (1999), ICEO)
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
Data =[0.806 9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020
0.517 9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374
0.100 8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982
0.908 2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426
0.965 8.074 8.777 3.467 1.863 6.708 6.349 4.534 0.276 7.633 1.567
0.669 7.650 5.658 0.720 2.764 3.278 5.283 7.474 6.274 1.409 8.208
0.524 1.256 3.605 8.623 6.905 4.584 8.133 6.071 6.888 4.187 5.448
0.902 8.314 2.261 4.224 1.781 4.124 0.932 8.129 8.658 1.208 5.762
0.531 0.226 8.858 1.420 0.945 1.622 4.698 6.228 9.096 0.972 7.637
0.876 7.305 2.228 1.242 5.928 9.133 1.826 4.060 5.204 8.713 8.247
0.462 0.652 7.027 0.508 4.876 8.807 4.632 5.808 6.937 3.291 7.016
0.491 2.699 3.516 5.874 4.119 4.461 7.496 8.817 0.690 6.593 9.789
0.463 8.327 3.897 2.017 9.570 9.825 1.150 1.395 3.885 6.354 0.109
0.714 2.132 7.006 7.136 2.641 1.882 5.943 7.273 7.691 2.880 0.564
0.352 4.707 5.579 4.080 0.581 9.698 8.542 8.077 8.515 9.231 4.670
0.869 8.304 7.559 8.567 0.322 7.128 8.392 1.472 8.524 2.277 7.826
0.813 8.632 4.409 4.832 5.768 7.050 6.715 1.711 4.323 4.405 4.591
0.811 4.887 9.112 0.170 8.967 9.693 9.867 7.508 7.770 8.382 6.740
0.828 2.440 6.686 4.299 1.007 7.008 1.427 9.398 8.480 9.950 1.675
0.964 6.306 8.583 6.084 1.138 4.350 3.134 7.853 6.061 7.457 2.258
0.789 0.652 2.343 1.370 0.821 1.310 1.063 0.689 8.819 8.833 9.070
0.360 5.558 1.272 5.756 9.857 2.279 2.764 1.284 1.677 1.244 1.234
0.369 3.352 7.549 9.817 9.437 8.687 4.167 2.570 6.540 0.228 0.027
0.992 8.798 0.880 2.370 0.168 1.701 3.680 1.231 2.390 2.499 0.064
0.332 1.460 8.057 1.336 7.217 7.914 3.615 9.981 9.198 5.292 1.224
0.817 0.432 8.645 8.774 0.249 8.081 7.461 4.416 0.652 4.002 4.644
0.632 0.679 2.800 5.523 3.049 2.968 7.225 6.730 4.199 9.614 9.229
0.883 4.263 1.074 7.286 5.599 8.291 5.200 9.214 8.272 4.398 4.506
0.608 9.496 4.830 3.150 8.270 5.079 1.231 5.731 9.494 1.883 9.732
0.326 4.138 2.562 2.532 9.661 5.611 5.500 6.886 2.341 9.699 6.500];
%
c   = Data(:,1);
dim = size(x,1);
A   = Data(:,2:dim+1);
f   = -sum(1./(c+sum((repmat(x',30,1)-A).^2,2)),1);
end
%
% End of shekel_foxholes.
