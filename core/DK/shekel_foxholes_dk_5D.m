function PB = shekel_foxholes_dk_5D()
%SHEKEL_FOXHOLES_DK_5D  DFO Benchmark Suite v2 core: shekel_foxholes, dimension 5, pattern DK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -10.40561723899244 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:K x3:D x4:K x5:D
%
lb_orig = [0 0 0 0 0];
ub_orig = [10 10 10 10 10];
lb_C_work = [0 0 0 0 0];
ub_C_work = [10 10 10 10 10];
T_break = [8.125 9.375 5 7.5 4.375];   % breakpoints of T (NaN for continuous)
T_image = [8.024914889028711 9.1517257530637934 5.1139278077814589 7.6208609576192838 4.5640883928378626];   % their images
m = 17;
ProblemData = cell(1, 5);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = {linspace(lb_orig(5), ub_orig(5), m)'};
PB.name        = 'shekel_foxholes_dk_5D';
PB.base_function = 'shekel_foxholes';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: DK adaptation of shekel_foxholes (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'D','K','D','K','D'};
PB.type_counts = struct('C', 0, 'D', 3, 'K', 2);
PB.proportion  = [0 3 2];
PB.pattern     = 'DK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [17 17 17 17 17];
PB.n_finite    = 5;
PB.finite_part_size = 1419857;
PB.M_K = 289;
PB.domain_size = 1419857;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shekel_foxholes_dk_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([2 13 17 15 10 5 11 4 3 12 8 14 7 6 9 16 1])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([4 6 5 11 13 16 10 9 12 15 1 2 17 8 14 3 7])';
D.fn = @(xh) shekel_foxholes_orig(xh);
PB.f_star      = -10.405617238992441;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [7 6 1 2 9];   % protocol starting point: Halton point 1050 (index representation)
PB.x0_alt_idx  = [11 5 12 3 4;3 6 5 1 17;16 6 16 2 12;8 5 8 3 8;12 6 3 3 3;4 6 11 1 15;14 5 7 2 10;6 5 14 3 7;10 6 10 1 2];   % 9 alternative starting points (Halton indices 1050 + 63 k)
PB.core_pidx   = 50;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shekel_foxholes_dk_5D: n = 5, pattern DK, m = 17, f* = -10.40561724 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('shekel_foxholes_dk_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.lb(1) + (k - 1) * (D.ub(1) - D.lb(1)) / (D.m - 1);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
k = check_idx(z(5), D.m, D.name, 5); g = D.lb(5) + (k - 1) * (D.ub(5) - D.lb(5)) / (D.m - 1);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
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
if ~iscell(xr) || numel(xr) ~= 5, error('shekel_foxholes_dk_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_foxholes_dk_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('shekel_foxholes_dk_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shekel_foxholes_dk_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shekel_foxholes_dk_5D: x%d unknown label %s.', i, lab); end
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
