function PB = neumaier3_ck_10D()
%NEUMAIER3_CK_10D  DFO Benchmark Suite v2.1.0 core: neumaier3, dimension 10, pattern CK, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -210 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C x10:K
%
lb_orig = [-100 -100 -100 -72.698042482299044 -66.211597568472087 -71.1692167641189 -100 -100 -100 -100];
ub_orig = [61.813232421875 69.889549357313427 72.739328 100 100 100 62.208259000730351 69.966651575724555 61.652338292101582 65.622042724178925];
lb_C_work = [-100 -100 -100 -72.698042482299044 -6621159756.847208 -71.1692167641189 -10000000000 -100 -10000000000 -100];
ub_C_work = [61.813232421875 69.889549357313427 72.739328 100 10000000000 100 6220825900.0730343 69.966651575724555 6165233829.2101574 65.622042724178925];
T_break = [NaN 16.799065183152976 NaN 24.44460641399418 NaN 25.113467665697996 NaN 27.474988681793405 NaN 13.865154372873008];   % breakpoints of T (NaN for continuous)
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
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CK adaptation of neumaier3 (F = f o T).';
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
D.cw = [10 NaN 24 NaN 3000000000.0000019 NaN 2799999999.9999981 NaN 1800000000 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [10 NaN 24 NaN 30 NaN 28 NaN 18 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([6 5 10 12 14 3 16 11 17 1 2 7 4 15 8 13 9])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([12 6 13 7 11 8 10 1 15 2 17 3 4 9 14 5 16])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([3 6 1 10 2 15 14 16 9 13 8 11 4 17 5 7 12])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([16 15 4 6 12 2 9 10 8 5 1 7 3 13 11 17 14])';
g = linspace(lb_orig(10), ub_orig(10), m); D.E{10} = g([14 15 5 17 10 6 9 13 8 7 12 4 3 2 16 1 11])';
D.fn = @(xh) neumaier3_orig(xh);
PB.f_star      = -210;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [31.552261710166931 2 -83.914513776640007 15 -4186053706.8207574 11 4827545108.3305511 9 -9326393231.5764313 12];   % protocol starting point: Halton point 1035 (index representation)
PB.x0_alt_idx  = [-48.09018862247467 2 36.312058511360021 15 7090360463.3017292 9 281223662.2893219 14 1981431632.4334011 15;-8.2689634561538696 1 -49.090265251839995 16 3381506633.2614422 7 -5215962399.9173975 3 -2143911594.3026485 2;-87.911413788795471 2 71.136307036160019 17 -1825873996.7951221 4 -9762283845.9586277 8 -6269254821.038702 5;41.981630206108093 3 -15.592654766080003 15 -5534727826.8354073 2 1856093182.813406 13 5740077683.4595833 7;-37.660820126533508 2 -94.361788334080003 16 5605456634.1946888 16 -2690228263.2278252 2 1614734456.7235298 10;2.1604050397872925 2 19.23159375872001 16 1896602804.1544046 13 -7236549709.2690563 8 -2510608770.0125208 13;-77.482045292854309 2 -59.537539809279998 15 -3310777825.9021592 11 3430962703.3374939 12 -6635951996.7485723 16;21.596955418586731 1 54.055842283520008 16 9464163144.2366066 7 -1115358742.7037373 1 5373380507.7497101 3;-58.045494914054871 2 -24.713291284479993 16 4256782514.1800404 5 -5661680188.7449694 6 546529640.52520943 5];   % 9 alternative starting points (Halton indices 1035 + 63 k)
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
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('neumaier3_ck_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('neumaier3_ck_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('neumaier3_ck_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('neumaier3_ck_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); tl = 1e-12 * (D.ubw(9) - D.lbw(9)); if v < D.lbw(9) - tl || v > D.ubw(9) + tl, error('neumaier3_ck_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
v = min(max(v, D.lbw(9)), D.ubw(9));
if abs(v - D.cw(9)) <= 1e-12 * (D.ubw(9) - D.lbw(9)), xh(9) = D.cx(9); else, xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9)); end
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
if ~iscell(xr) || numel(xr) ~= 10, error('neumaier3_ck_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('neumaier3_ck_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('neumaier3_ck_10D: x%d must be a grid value.', i); end
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

