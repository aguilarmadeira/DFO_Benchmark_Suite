function PB = paviani_c_10D()
%PAVIANI_C_10D  DFO Benchmark Suite v2 core: paviani, dimension 10, pattern C, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -45.77845205382886 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
%
lb_orig = [2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999 2.0009999999999999];
ub_orig = [9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006];
lb_C_work = [2.0009999999999994 2.0009999999999994 2.0009999999999994 2.0009999999999994 2.0009999999999994 200099999.99999994 200099999.99999994 200099999.99999994 200099999.99999994 200099999.99999994];
ub_C_work = [9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 9.9990000000000006 999900000 999900000 999900000 999900000 999900000];
T_break = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 10);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = [lb_C_work(6); ub_C_work(6)];
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = [lb_C_work(8); ub_C_work(8)];
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
ProblemData{10} = [lb_C_work(10); ub_C_work(10)];
PB.name        = 'paviani_c_10D';
PB.base_function = 'paviani';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of paviani (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 10, 'D', 0, 'K', 0);
PB.proportion  = [10 0 0];
PB.pattern     = 'C';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 1 1 1 1 100000000 100000000 100000000 100000000 100000000];
PB.cardinality = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'paviani_c_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.fn = @(xh) paviani_orig(xh);
PB.f_star      = -45.778452053828865;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 1;
PB.x0_idx      = [3.5045302734374997 5.1753319615912208 4.3453737600000002 2.1742178259058722 3.8998489857250185 760724487.93809724 978248585.38571143 656728779.70549631 268201652.00953394 785004374.92312098];   % protocol starting point: Halton point 1036 (index representation)
PB.x0_alt_idx  = [8.5657646484375007 5.4057270233196153 9.9119817600000015 2.5239845897542685 9.3259902329075892 661341055.98543465 754083187.46183586 873965607.23137474 827673814.41604328 924802948.05035067;4.5355224609375 4.8461961591220852 5.9577705600000002 2.8737513536026649 7.5413125469571742 561957624.03277183 483033523.30551589 335596078.14550221 623566063.94345355 265752531.87912577;6.5975068359374998 4.9888216735253774 2.3106825599999996 2.1042644731361926 5.0355529676934632 401415157.0323168 258868125.38164049 594811036.59425569 419458313.47086376 405551105.00635523;2.5204013671874996 5.4825253772290807 7.5087427200000008 2.4540312369845894 3.25087528174305 302031725.07965404 831735253.40932226 854025995.04300904 215350562.99827397 545349678.13358474;9.0812607421874993 4.9229945130315498 3.8616547199999998 2.8037980008329857 8.6114636295335032 997715748.74829316 607569855.48544669 315656465.95713651 809530903.26292419 685148251.26081419;5.0510185546875004 5.1533895747599452 9.1211395199999998 2.9936713869221152 6.8267859435830882 837173281.74783802 383404457.5615713 574871424.40588999 605423152.79033446 824946824.38804364;7.1130029296875001 5.2667585733882021 5.4740515199999997 2.2241845064556429 4.3210263643193763 737789849.79517531 909387319.35680842 792108251.93176842 401315402.31774467 964745397.51527309;3.0124658203124994 4.707227709190672 3.0554563199999998 2.5739512703040397 2.5363486783689635 577247382.79471993 685221921.43293309 253738722.84589583 995495742.58239484 305694981.3440485;8.0737001953125009 5.2009314128943753 7.0864483200000006 2.923718034152436 7.9624899255515329 477863950.84205729 461056523.50905752 512953681.2946493 756679814.25166416 445493554.47127795];   % 9 alternative starting points (Halton indices 1036 + 63 k)
PB.core_pidx   = 36;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('paviani_c_10D: n = 10, pattern C, m = 0, f* = -45.77845205 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('paviani_c_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('paviani_c_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('paviani_c_10D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('paviani_c_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('paviani_c_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('paviani_c_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
v = z(6); if v < D.lbw(6) || v > D.ubw(6), error('paviani_c_10D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('paviani_c_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
v = z(8); if v < D.lbw(8) || v > D.ubw(8), error('paviani_c_10D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('paviani_c_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
v = z(10); if v < D.lbw(10) || v > D.ubw(10), error('paviani_c_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
xh(10) = map_work_to_orig(v, D.lbw(10), D.ubw(10), D.lb(10), D.ub(10));
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
if ~iscell(xr) || numel(xr) ~= 10, error('paviani_c_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('paviani_c_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('paviani_c_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('paviani_c_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('paviani_c_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = paviani_orig(x)
%
% Purpose:
%
%    Function paviani is the function described in
%    Montaz et al. (2005).
%
%    dim = 10
%    Minimum global value = -45.778
%    Global minima = 9.351*ones(10,1)
%    Local minima
%    Search domain: 2.001 <= x <= 9.999
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
f = sum((log(x-2)).^2+(log(10-x)).^2,1)-prod(x,1)^0.2;
end
%
% End of paviani.

