function PB = shekel_foxholes_cdk_5D()
%SHEKEL_FOXHOLES_CDK_5D  DFO Benchmark Suite v2.1.0 core: shekel_foxholes, dimension 5, pattern CDK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -10.40561723899244 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C x5:D
%
lb_orig = [0 1.3797337804196517 0 1.4936693044564764 1.9497507000887921];
ub_orig = [8.4762084960937507 10 8.6873664000000002 10 10];
lb_C_work = [0 1.3797337804196517 0 1.4936693044564764 1.9497507000887921];
ub_C_work = [8.4762084960937507 10 8.6873664000000002 10 10];
T_break = [NaN 9.1379733780419645 5.2124198399999999 NaN 4.5660817225599342];   % breakpoints of T (NaN for continuous)
T_image = [NaN 9.1517257530637934 5.1139278077814589 NaN 4.5640883928378626];   % their images
m = 41;
ProblemData = cell(1, 5);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = {linspace(lb_orig(5), ub_orig(5), m)'};
PB.name        = 'shekel_foxholes_cdk_5D';
PB.base_function = 'shekel_foxholes';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CDK adaptation of shekel_foxholes (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C','D'};
PB.type_counts = struct('C', 2, 'D', 2, 'K', 1);
PB.proportion  = [2 2 1];
PB.pattern     = 'CDK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1];
PB.cardinality = [Inf 41 41 Inf 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 41;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shekel_foxholes_cdk_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
D.cw = [8.024914889028711 NaN NaN 7.6208609576192838 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [8.024914889028711 NaN NaN 7.6208609576192838 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([7 10 31 29 12 3 22 11 1 37 30 5 23 4 36 13 38 28 18 20 6 27 26 34 33 8 41 17 40 35 32 15 9 21 25 16 19 14 24 39 2])';
D.fn = @(xh) shekel_foxholes_orig(xh);
PB.f_star      = -10.405617238992441;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [2.9178354442119603 14 2 2.0250935461605701 22];   % protocol starting point: Halton point 1050 (index representation)
PB.x0_alt_idx  = [5.1031079471111305 11 30 2.3970905153534354 9;0.83189350962638864 12 10 1.5786971831291314 40;7.7850332915782934 13 38 1.9506941523219967 27;3.4641535699367525 10 18 2.322691121514862 18;5.6494260728359231 13 7 2.6946880907077277 6;1.3782116353511811 14 26 1.7062390011381139 36;6.7420623242855076 10 15 2.0782359703309794 24;2.3963499605655674 12 34 2.4502329395238447 15;4.5816224634647371 14 23 1.6318396072995407 2];   % 9 alternative starting points (Halton indices 1050 + 63 k)
PB.core_pidx   = 50;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shekel_foxholes_cdk_5D: n = 5, pattern CDK, m = 41, f* = -10.40561724 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('shekel_foxholes_cdk_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('shekel_foxholes_cdk_5D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('shekel_foxholes_cdk_5D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
k = check_idx(z(5), D.m, D.name, 5); g = D.lb(5) + (k - 1) * (D.ub(5) - D.lb(5)) / (D.m - 1);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
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
if ~iscell(xr) || numel(xr) ~= 5, error('shekel_foxholes_cdk_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_foxholes_cdk_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_foxholes_cdk_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shekel_foxholes_cdk_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shekel_foxholes_cdk_5D: x%d unknown label %s.', i, lab); end
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
