function PB = neumaier3_ck_10D()
%NEUMAIER3_CK_10D  DFO Benchmark Suite v2 core: neumaier3, dimension 10, pattern CK, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K
%
lb_orig = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub_orig = [100 100 100 100 100 100 100 100 100 100];
lb_C_work = [-100 -100 -100 -100 -10000000000 -100 -10000000000 -100 -10000000000 -100];
ub_C_work = [100 100 100 100 10000000000 100 10000000000 100 10000000000 100];
T_break = [NaN 12.5 NaN 25 NaN 25 NaN 25 NaN 12.5];   % breakpoints of T (NaN for continuous)
T_image = [NaN 18 NaN 28 NaN 30 NaN 24 NaN 10];   % their images
m = 17;
ProblemData = cell(1, 10);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
ProblemData{10} = arrayfun(@(j) sprintf('k10_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'neumaier3_ck_10D';
PB.base_function = 'neumaier3';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of neumaier3 (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C','K','C','K','C','K'};
PB.type_counts = struct('C', 5, 'D', 0, 'K', 5);
PB.proportion  = [5 0 5];
PB.pattern     = 'CK';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 1 100000000 100000000 100000000];
PB.cardinality = [Inf 17 Inf 17 Inf 17 Inf 17 Inf 17];
PB.n_finite    = 5;
PB.finite_part_size = 1419857;
PB.M_K = 1419857;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'neumaier3_ck_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([6 5 10 12 14 3 16 11 17 1 2 7 4 15 8 13 9])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([12 6 13 7 11 8 10 1 15 2 17 3 4 9 14 5 16])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([3 6 1 10 2 15 14 16 9 13 8 11 4 17 5 7 12])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([16 15 4 6 12 2 9 10 8 5 1 7 3 13 11 17 14])';
g = linspace(lb_orig(10), ub_orig(10), m); D.E{10} = g([14 15 5 17 10 6 9 13 8 7 12 4 3 2 16 1 11])';
D.fn = @(xh) neumaier3_orig(xh);
PB.f_star      = -210;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [62.59765625 2 -81.376000000000005 15 -7069872276.4838467 11 8282108691.227356 9 -9166598175.3924541 12];   % protocol starting point: Halton point 1035 (index representation)
PB.x0_alt_idx  = [-35.83984375 2 57.824000000000041 15 6498873027.7986488 9 2676572359.0474262 14 4823703460.1791725 15;13.37890625 1 -41.05599999999999 16 2036063110.4432755 7 -4101363728.8825569 3 -280266294.07413483 2;-85.05859375 2 98.144000000000005 17 -4229902329.0758839 4 -9706900061.0624866 8 -5384236048.3274431 5;75.48828125 3 -2.2719999999999914 15 -8692712246.4312553 2 4618359454.5084476 13 9473987014.054409 7;-22.94921875 2 -93.471999999999994 16 4712109828.5636215 16 -987176877.67148399 2 4370017259.8010998 10;26.26953125 2 38.048000000000002 16 249299911.20825195 13 -6592713209.8514147 8 -733952494.45220757 13;-72.16796875 2 -53.152000000000001 15 -6016665528.3109074 11 6560146549.9694691 12 -5837922248.7055149 16;50.29296875 1 78.367999999999995 16 9355235298.135376 7 954610217.78953934 1 9020300813.6763382 3;-48.14453125 2 -12.831999999999994 16 3089269858.6162148 5 -4650926114.3903933 6 3048409632.6128044 5];   % 9 alternative starting points (Halton indices 1035 + 63 k)
PB.core_pidx   = 35;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('neumaier3_ck_10D: n = 10, pattern CK, m = 17, f* = -210 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('neumaier3_ck_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('neumaier3_ck_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('neumaier3_ck_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('neumaier3_ck_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('neumaier3_ck_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('neumaier3_ck_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
k = check_idx(z(10), D.m, D.name, 10); g = D.E{10}(k);
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
if ~iscell(xr) || numel(xr) ~= 10, error('neumaier3_ck_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('neumaier3_ck_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('neumaier3_ck_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('neumaier3_ck_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('neumaier3_ck_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = neumaier3_orig(x)
%
% Purpose:
%
%    Function neumaier3 is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -n*(n+4)*(n-1)/6
%    Global minima = [1:n]'.*(n+1-[1:n]')
%    Local minima (several)
%    Search domain: -n^2 <= x <= n^2   (Montaz et al. (2005))
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
dim  = size(x,1);
xaux = x(1:dim-1,1);
f    = sum((x-1).^2,1)-sum(x(2:dim,1).*xaux,1);
end
%
% End of neumaier3.

