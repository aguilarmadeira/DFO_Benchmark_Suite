function PB = griewank_cd_10D()
%GRIEWANK_CD_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern CD, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C x10:D
%
lb_orig = [-270.0029296875 -256.74338261443881 -400 -400 -400 -400 -400 -400 -278.2615270814498 -287.59161917257779];
ub_orig = [400 400 264.52531199999999 285.07788421491045 275.49101837306193 285.00822800322118 291.6565654146861 262.337520430322 400 400];
lb_C_work = [-270.0029296875 -256.74338261443881 -400 -400 -40000000000 -400 -40000000000 -400 -27826152708.144981 -287.59161917257779];
ub_C_work = [400 400 264.52531199999999 285.07788421491045 27549101837.30619 285.00822800322118 29165656541.468613 262.337520430322 40000000000 400];
T_break = [NaN -10.464614134024259 NaN -14.643690129112883 NaN -14.68287174818812 NaN 13.960950268951251 NaN 13.22971421542502];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 NaN 0 NaN 0 NaN 0 NaN 0];   % their images
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
PB.name        = 'griewank_cd_10D';
PB.base_function = 'griewank';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CD adaptation of griewank (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','D','C','D','C','D','C','D','C','D'};
PB.type_counts = struct('C', 5, 'D', 5, 'K', 0);
PB.proportion  = [5 5 0];
PB.pattern     = 'CD';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 1 100000000 100000000 100000000];
PB.cardinality = [Inf 17 Inf 17 Inf 17 Inf 17 Inf 17];
PB.n_finite    = 5;
PB.finite_part_size = 1419857;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'griewank_cd_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.cw = [0 NaN 5.6843418860808015e-14 NaN 0 NaN 0 NaN 0 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [0 NaN 0 NaN 0 NaN 0 NaN 0 NaN];
D.fn = @(xh) griewank_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [64.344235420227051 15 -24.463455682559982 3 24047307677.445847 11 -29666681884.502754 14 4361585085.6274414 5];   % protocol starting point: Halton point 1022 (index representation)
PB.x0_alt_idx  = [222.35762166976929 14 -327.48699795456002 1 8974367598.0469437 9 19874117091.566452 3 -15890935074.534527 8;-107.40944528579712 14 109.50484721664003 2 -12188649281.1091 7 488587057.45241547 8 34497870245.505127 11;308.72518682479858 15 -193.51869505535998 3 -33351666160.265129 4 -22951432918.437107 13 17188738647.018219 13;-21.041880130767822 14 238.36959571967998 1 18566238557.664421 2 26589366057.632111 2 -120392951.46868134 16;139.91585493087769 14 -64.653946552320008 2 -3150421666.9240875 16 7203836023.5180664 7 -17429524549.955582 3;-189.85121202468872 15 -393.19526080511997 2 -18223361746.322987 12 -12181694010.595978 13 32959280770.084068 6;351.90896940231323 13 69.314356346880061 3 -39386378625.479019 10 -35621713986.485497 1 15650149171.597168 9;22.141902446746826 14 -259.22695790592002 1 12531526092.450539 7 13919084989.58371 6 -1658982426.8897362 11;181.13673830032349 15 203.28265924608002 2 -8631490786.7055016 5 -5466445044.5303268 11 -21911502587.051712 14];   % 9 alternative starting points (Halton indices 1022 + 63 k)
PB.core_pidx   = 22;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('griewank_cd_10D: n = 10, pattern CD, m = 17, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('griewank_cd_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('griewank_cd_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('griewank_cd_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
k = check_idx(z(4), D.m, D.name, 4); g = D.lb(4) + (k - 1) * (D.ub(4) - D.lb(4)) / (D.m - 1);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('griewank_cd_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
k = check_idx(z(6), D.m, D.name, 6); g = D.lb(6) + (k - 1) * (D.ub(6) - D.lb(6)) / (D.m - 1);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('griewank_cd_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
k = check_idx(z(8), D.m, D.name, 8); g = D.lb(8) + (k - 1) * (D.ub(8) - D.lb(8)) / (D.m - 1);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); tl = 1e-12 * (D.ubw(9) - D.lbw(9)); if v < D.lbw(9) - tl || v > D.ubw(9) + tl, error('griewank_cd_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
v = min(max(v, D.lbw(9)), D.ubw(9));
if abs(v - D.cw(9)) <= 1e-12 * (D.ubw(9) - D.lbw(9)), xh(9) = D.cx(9); else, xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9)); end
k = check_idx(z(10), D.m, D.name, 10); g = D.lb(10) + (k - 1) * (D.ub(10) - D.lb(10)) / (D.m - 1);
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
if ~iscell(xr) || numel(xr) ~= 10, error('griewank_cd_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_cd_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_cd_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('griewank_cd_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('griewank_cd_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = griewank_orig(x)
%
% Purpose:
%
%    Function griewank is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -600 <= x <= 600 (Huyer and Neumaier (1999))
%                   -400 <= x <= 400 (Storn and Price (1997))
%    Cases considered: n = 10 Storn and Price (1997)
%                      n = 5 Huyer and Neumaier (1999), ICEO
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
f = 1 + 1/4000*sum(x.^2,1)- prod(cos(x./sqrt([1:n]')),1);
end
%
% End of griewank.

