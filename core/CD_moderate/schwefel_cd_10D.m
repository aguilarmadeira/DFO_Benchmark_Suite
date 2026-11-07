function PB = schwefel_cd_10D()
%SCHWEFEL_CD_10D  DFO Benchmark Suite v2 core: schwefel, dimension 10, pattern CD, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -4189.828872721624 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D
%
lb_orig = [-500 -500 -500 -500 -500 -500 -500 -500 -500 -500];
ub_orig = [500 500 500 500 500 500 500 500 500 500];
lb_C_work = [-500 -500 -50000 -500 -5 -500 -500 -500 -50000 -500];
ub_C_work = [500 500 50000 500 5 500 500 500 50000 500];
T_break = [NaN 437.5 NaN 437.5 NaN 437.5 NaN 437.5 NaN 437.5];   % breakpoints of T (NaN for continuous)
T_image = [NaN 420.96870000000001 NaN 420.96870000000001 NaN 420.96870000000001 NaN 420.96870000000001 NaN 420.96870000000001];   % their images
m = 17;
ProblemData = cell(1, 10);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = {linspace(lb_orig(2), ub_orig(2), m)'};
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = {linspace(lb_orig(4), ub_orig(4), m)'};
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = {linspace(lb_orig(6), ub_orig(6), m)'};
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = {linspace(lb_orig(8), ub_orig(8), m)'};
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
ProblemData{10} = {linspace(lb_orig(10), ub_orig(10), m)'};
PB.name        = 'schwefel_cd_10D';
PB.base_function = 'schwefel';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CD adaptation of schwefel (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','D','C','D','C','D','C','D','C','D'};
PB.type_counts = struct('C', 5, 'D', 5, 'K', 0);
PB.proportion  = [5 5 0];
PB.pattern     = 'CD';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01 1 100];
PB.cardinality = [Inf 17 Inf 17 Inf 17 Inf 17 Inf 17];
PB.n_finite    = 5;
PB.finite_part_size = 1419857;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'schwefel_cd_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.fn = @(xh) schwefel_orig(xh);
PB.f_star      = -4189.8288727216241;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-93.26171875 13 -12688 9 -3.45229151014275 9 64.624465703236297 3 1993.0960795594656 2];   % protocol starting point: Halton point 1046 (index representation)
PB.x0_alt_idx  = [164.55078125 12 37712.000000000015 9 3.332081141998497 7 -215.65235090576022 7 -27866.359825758198 5;-339.35546875 13 7472 10 1.100676183320811 5 441.99063708528399 12 46424.755486151058 8;293.45703125 13 -38128 8 -2.0323065364387682 1 161.71382047628742 1 20904.906714884521 11;-216.30859375 12 26864.000000000015 9 -4.2637114951164534 15 -118.5629961327092 7 -4614.9420563820167 13;41.50390625 13 -18735.999999999996 10 2.4386995423809843 12 -398.83981274170571 12 -30134.790827648554 16;-462.40234375 13 47024.000000000015 8 0.20729458370329912 10 317.42316303684106 17 44156.324484260695 3;451.66015625 12 1424 9 -2.9256881360562801 8 -21.473641359658075 5 18636.475712994157 6;-61.03515625 12 -48016 9 4.760262277166861 5 -301.75045796865459 10 -6883.3730582723801 9;196.77734375 13 21584 10 1.6272795574072809 3 414.51251780989219 16 -36742.828963590044 11];   % 9 alternative starting points (Halton indices 1046 + 63 k)
PB.core_pidx   = 46;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('schwefel_cd_10D: n = 10, pattern CD, m = 17, f* = -4189.828873 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('schwefel_cd_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('schwefel_cd_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('schwefel_cd_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.lb(4) + (k - 1) * (D.ub(4) - D.lb(4)) / (D.m - 1);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('schwefel_cd_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.lb(6) + (k - 1) * (D.ub(6) - D.lb(6)) / (D.m - 1);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('schwefel_cd_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.lb(8) + (k - 1) * (D.ub(8) - D.lb(8)) / (D.m - 1);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('schwefel_cd_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
k = check_idx(z(10), D.m, D.name, 10); g = D.lb(10) + (k - 1) * (D.ub(10) - D.lb(10)) / (D.m - 1);
xh(10) = T_map(g, D.lb(10), D.ub(10), D.Tb(10), D.Ti(10));
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
if ~iscell(xr) || numel(xr) ~= 10, error('schwefel_cd_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('schwefel_cd_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('schwefel_cd_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('schwefel_cd_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('schwefel_cd_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = schwefel_orig(x)
%
% Purpose:
%
%    Function schwefel is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -418.9829*n
%    Global minima = 420.97*ones(n,1)
%    Local minima 
%    Search domain: -500 <= x <= 500 (Montaz et al. (2005))
%    Cases considered: n = 10 Montaz et al. (2005)
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
f = -sum(x.*sin(sqrt(abs(x))),1);
end
%
% End of schwefel.

