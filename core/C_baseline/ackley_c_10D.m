function PB = ackley_c_10D()
%ACKLEY_C_10D  DFO Benchmark Suite v2 core: ackley, dimension 10, pattern C, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 4.440892098500626e-16 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
%
lb_orig = [-30 -30 -30 -30 -30 -30 -30 -30 -30 -30];
ub_orig = [30 30 30 30 30 30 30 30 30 30];
lb_C_work = [-30 -30 -30 -30 -30 -30 -30 -30 -30 -30];
ub_C_work = [30 30 30 30 30 30 30 30 30 30];
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
PB.name        = 'ackley_c_10D';
PB.base_function = 'ackley';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of ackley (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 10, 'D', 0, 'K', 0);
PB.proportion  = [10 0 0];
PB.pattern     = 'C';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1 1 1 1 1 1 1 1];
PB.cardinality = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'ackley_c_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.fn = @(xh) ackley_orig(xh);
PB.f_star      = 4.4408920985006262e-16;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [5.56640625 10.850480109739365 -17.692799999999998 -25.22698875468555 -28.151765589782119 -25.603095129722348 24.431101160187261 13.396996646741499 3.5777101997205563 1.3936610767149169];   % protocol starting point: Halton point 1001 (index representation)
PB.x0_alt_idx  = [-25.283203125 12.578875171467764 24.0672 -23.802582257392753 17.963936889556727 21.998179335457444 7.6144921636474692 -26.990814987607521 -14.337963343470042 11.881175939972941;23.935546875 16.28257887517146 -5.5967999999999982 -29.575177009579342 -0.83395942900074971 14.54255803368229 -9.202116832892326 -10.693978714098268 -29.649872606229966 22.368690803230969;-6.298828125 12.085048010973928 -21.436800000000002 -26.951270304039983 -14.222389181066866 2.4988620846609066 -29.53592509668227 8.7520046654031205 14.924796580915583 -27.072450694985445;9.169921875 13.155006858710571 6.0384000000000029 -24.327363598500625 26.483846731780623 -4.9567592171142465 13.439853450030533 28.197988044904506 -0.38711268184433223 -16.584935831727421;-21.416015625 14.883401920438956 -9.8016000000000005 -21.703456892961267 13.0954169797145 -12.412380518889393 -3.3767555465092656 -12.189823589444526 -15.699021944604256 -6.0974209684693932;16.552734375 10.685871056241425 18.134400000000007 -28.67555185339442 -6.1942490267058261 -24.456076467910787 -20.193364543049054 7.25615979005687 28.875647242541298 4.3900938947886345;-13.681640625 14.389574759945134 2.2944000000000031 -26.051645147855062 -19.582678778771943 27.733272644515246 22.782414003663746 26.702143169558241 13.563737979781372 14.877608758046662;1.787109375 15.240054869684499 -27.369599999999998 -23.4277384423157 21.123557134075547 20.2776513427401 2.4486057398738126 -16.834815570782915 -1.7481712829785465 25.36512362130469;-28.974609375 11.042524005486968 14.3904 -29.200333194502292 7.7351273820094235 8.2339553937187091 -14.36800325666599 2.6111678087184629 -17.060080545738472 -24.076017876911724];   % 9 alternative starting points (Halton indices 1001 + 63 k)
PB.core_pidx   = 1;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('ackley_c_10D: n = 10, pattern C, m = 0, f* = 4.440892099e-16 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('ackley_c_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('ackley_c_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('ackley_c_10D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('ackley_c_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('ackley_c_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('ackley_c_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
v = z(6); if v < D.lbw(6) || v > D.ubw(6), error('ackley_c_10D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('ackley_c_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
v = z(8); if v < D.lbw(8) || v > D.ubw(8), error('ackley_c_10D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('ackley_c_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
v = z(10); if v < D.lbw(10) || v > D.ubw(10), error('ackley_c_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
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
if ~iscell(xr) || numel(xr) ~= 10, error('ackley_c_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('ackley_c_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('ackley_c_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('ackley_c_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('ackley_c_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = ackley_orig(x)
%
% Purpose:
%
%    Function ackley is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima
%    Search domain: -30 <= x <= 30 (Montaz et al. (2005))
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
n = size(x,1);
f = -20*exp(-0.02*sqrt(1/n*sum(x.^2,1)))-exp(1/n*sum(cos(2*pi.*x),1))+20+exp(1);
end
%
% End of ackley.

