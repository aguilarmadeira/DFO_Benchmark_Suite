function PB = shekel_foxholes_k_10D()
%SHEKEL_FOXHOLES_K_10D  DFO Benchmark Suite v2.1.0 core: shekel_foxholes, dimension 10, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 9 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -10.20879279215383 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K x3:K x4:K x5:K x6:K x7:K x8:K x9:K x10:K
%
lb_orig = [0 1.6130671137529848 1.4526336 1.5936693044564765 1.3191721876920974 1.6072196351668362 1.8508327247039666 0 0 0];
ub_orig = [8.1262084960937493 10 10 10 10 10 10 8.6069198364039572 8.2943371414481799 8.594062897207758];
lb_C_work = [0 1.6130671137529848 1.4526336 1.5936693044564765 1.3191721876920974 1.6072196351668362 1.8508327247039666 0 0 0];
ub_C_work = [8.1262084960937493 10 10 10 10 10 10 8.6069198364039572 8.2943371414481799 8.594062897207758];
T_break = [8.1262084960937493 8.951633389219122 4.6578959999999991 7.8984173261141191 4.5744826173075612 4.754512271979273 2.8694786341159708 6.4551898773029679 1.0367921426810225 5.3712893107548485];   % breakpoints of T (NaN for continuous)
T_image = [8.0249668149171214 9.1519283831871547 5.1139910122822458 7.6209592438519174 4.5640220516540477 4.7110050663986964 2.9960306020294212 6.1259930245105814 0.73405773897170223 4.9819993796057354];   % their images
m = 9;
ProblemData = cell(1, 10);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = arrayfun(@(j) sprintf('k3_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = arrayfun(@(j) sprintf('k5_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = arrayfun(@(j) sprintf('k7_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{9} = arrayfun(@(j) sprintf('k9_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{10} = arrayfun(@(j) sprintf('k10_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'shekel_foxholes_k_10D';
PB.base_function = 'shekel_foxholes';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of shekel_foxholes (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'K','K','K','K','K','K','K','K','K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 10);
PB.proportion  = [0 0 10];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [9 9 9 9 9 9 9 9 9 9];
PB.n_finite    = 10;
PB.finite_part_size = 3486784401;
PB.M_K = 3486784401;
PB.domain_size = 3486784401;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shekel_foxholes_k_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.cw = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([7 6 3 1 2 9 8 5 4])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([1 3 5 4 9 8 7 2 6])';
g = linspace(lb_orig(3), ub_orig(3), m); D.E{3} = g([6 4 7 9 2 5 1 8 3])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([1 5 4 6 2 7 9 3 8])';
g = linspace(lb_orig(5), ub_orig(5), m); D.E{5} = g([7 1 4 5 8 6 9 2 3])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([4 6 1 8 2 3 5 9 7])';
g = linspace(lb_orig(7), ub_orig(7), m); D.E{7} = g([4 8 3 6 2 9 7 5 1])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([3 7 6 4 2 5 9 1 8])';
g = linspace(lb_orig(9), ub_orig(9), m); D.E{9} = g([8 9 7 4 2 1 3 6 5])';
g = linspace(lb_orig(10), ub_orig(10), m); D.E{10} = g([7 1 9 3 2 5 8 4 6])';
D.fn = @(xh) shekel_foxholes_orig(xh);
PB.f_star      = -10.208792792153835;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [8 6 3 3 6 8 8 4 7 3];   % protocol starting point: Halton point 1051 (index representation)
PB.x0_alt_idx  = [4 6 8 3 3 7 6 6 5 4;6 6 4 2 2 6 3 9 2 6;2 6 2 3 7 4 9 3 8 7;8 6 6 3 5 3 6 6 6 9;4 6 4 3 3 2 4 8 4 2;6 6 8 2 9 9 2 3 2 3;2 6 5 3 6 8 7 5 8 5;7 6 1 3 5 6 5 8 6 6;3 6 7 2 2 5 3 2 4 8];   % 9 alternative starting points (Halton indices 1051 + 63 k)
PB.core_pidx   = 51;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shekel_foxholes_k_10D: n = 10, pattern K, m = 9, f* = -10.20879279 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('shekel_foxholes_k_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.E{3}(k);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
k = check_idx(z(5), D.m, D.name, 5); g = D.E{5}(k);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
k = check_idx(z(7), D.m, D.name, 7); g = D.E{7}(k);
xh(7) = T_map(g, D.lb(7), D.ub(7), D.Tb(7), D.Ti(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
k = check_idx(z(9), D.m, D.name, 9); g = D.E{9}(k);
xh(9) = T_map(g, D.lb(9), D.ub(9), D.Tb(9), D.Ti(9));
k = check_idx(z(10), D.m, D.name, 10); g = D.E{10}(k);
xh(10) = T_map(g, D.lb(10), D.ub(10), D.Tb(10), D.Ti(10));
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
xr = cell(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 10, error('shekel_foxholes_k_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_foxholes_k_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_foxholes_k_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shekel_foxholes_k_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shekel_foxholes_k_10D: x%d unknown label %s.', i, lab); end
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
