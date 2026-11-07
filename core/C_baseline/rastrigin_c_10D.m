function PB = rastrigin_c_10D()
%RASTRIGIN_C_10D  DFO Benchmark Suite v2.1.0 core: rastrigin, dimension 10, pattern C, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
%
lb_orig = [-5.1200000000000001 -5.1200000000000001 -3.6955815936 -3.2101397750937108 -3.775091646745441 -3.3681715906305802 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001];
ub_orig = [3.7836375000000002 3.2332190011685213 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 3.6885590546090206 3.3936609817297292 3.7162186241472837 3.2362623805814099];
lb_C_work = [-5.1200000000000001 -5.1200000000000001 -3.6955815936 -3.2101397750937108 -3.775091646745441 -3.3681715906305802 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001 -5.1200000000000001];
ub_C_work = [3.7836375000000002 3.2332190011685213 5.1200000000000001 5.1200000000000001 5.1200000000000001 5.1200000000000001 3.6885590546090206 3.3936609817297292 3.7162186241472837 3.2362623805814099];
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
PB.name        = 'rastrigin_c_10D';
PB.base_function = 'rastrigin';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: C adaptation of rastrigin (F = f o T).';
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
D = struct('name', 'rastrigin_c_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.cw = [0 0 0 0 0 0 0 0 0 0];   % anchors of the continuous maps (NaN for finite)
D.cx = [0 0 0 0 0 0 0 0 0 0];
D.fn = @(xh) rastrigin_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-4.8374138488769534 1.9078568642661535 -2.522051371859968 1.7303512868527466 1.5713195609152306 -3.244538686260904 -3.2553752459203378 1.5330431348697102 -2.8308768715942936 2.1436337065203936];   % protocol starting point: Halton point 1040 (index representation)
PB.x0_alt_idx  = [3.3010673034667972 2.1484845721187442 1.9210017513144342 2.0946431137518884 -1.2154972810778952 4.1386638215656957 3.0538694748549817 -4.644605313310616 3.3501913218895192 -4.742085472722076;-1.1855312805175782 1.5641029959053094 -0.74483012259020764 2.4589349406510301 -3.2003524419219187 3.0839206061618949 0.068676959910138002 -1.8853301474868536 1.0952016919086418 -3.2814783740948856;1.1099377624511719 1.7130630055283422 3.6982230005841927 1.6574929214729179 2.8344092087250639 1.3801046428172965 -2.4001579308280312 0.87394501833690885 -1.1597879380722333 -1.8208712754676957;-3.4288305725097654 2.2286938080696093 0.96468746004070471 2.0217847483720606 0.84955404788104083 0.32536142741349572 -4.8689928215661995 -4.8568572491432134 -3.4147775680531089 -0.36026417684050482;2.2055025329589846 1.6443122318561736 -3.0552177466408956 2.3860765752712023 -2.0101684014892758 -1.378454535931104 1.4402518992091196 -2.097582083319451 3.1497477992245502 1.1003429217866865;-2.2810960510253908 1.8849399397087643 2.7419087093104659 2.5838349955878783 4.8265549303068074 -2.4331977513349039 -1.0285829915290492 0.66169308250431147 0.89475816924367546 2.5609500204138769;0.014372991943359281 2.003344049921945 -1.2779964973711357 1.7823929764097661 2.0397380883136833 4.9500047564916967 -4.0137755064738929 2.9741220676278681 -1.3602314607372006 -4.3247691588285928;-4.550480218505859 1.4189624737085085 4.5191299585802245 2.146684803308907 0.05488292746965806 3.2461887931470965 2.295469214301427 -2.7566801998522528 -3.6152210907180757 -2.8641620602014028;2.7532849182128905 1.9345932762497755 0.49922475189862414 2.5109766302080496 -2.7319339145234665 2.1914455777432957 -0.17336567643674172 0.0025949659715100992 2.5658471027657352 -1.403554961574212];   % 9 alternative starting points (Halton indices 1040 + 63 k)
PB.core_pidx   = 40;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('rastrigin_c_10D: n = 10, pattern C, m = 0, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('rastrigin_c_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('rastrigin_c_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
v = z(2); tl = 1e-12 * (D.ubw(2) - D.lbw(2)); if v < D.lbw(2) - tl || v > D.ubw(2) + tl, error('rastrigin_c_10D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
v = min(max(v, D.lbw(2)), D.ubw(2));
if abs(v - D.cw(2)) <= 1e-12 * (D.ubw(2) - D.lbw(2)), xh(2) = D.cx(2); else, xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2)); end
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('rastrigin_c_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('rastrigin_c_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('rastrigin_c_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
v = z(6); tl = 1e-12 * (D.ubw(6) - D.lbw(6)); if v < D.lbw(6) - tl || v > D.ubw(6) + tl, error('rastrigin_c_10D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
v = min(max(v, D.lbw(6)), D.ubw(6));
if abs(v - D.cw(6)) <= 1e-12 * (D.ubw(6) - D.lbw(6)), xh(6) = D.cx(6); else, xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6)); end
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('rastrigin_c_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
v = z(8); tl = 1e-12 * (D.ubw(8) - D.lbw(8)); if v < D.lbw(8) - tl || v > D.ubw(8) + tl, error('rastrigin_c_10D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
v = min(max(v, D.lbw(8)), D.ubw(8));
if abs(v - D.cw(8)) <= 1e-12 * (D.ubw(8) - D.lbw(8)), xh(8) = D.cx(8); else, xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8)); end
v = z(9); tl = 1e-12 * (D.ubw(9) - D.lbw(9)); if v < D.lbw(9) - tl || v > D.ubw(9) + tl, error('rastrigin_c_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
v = min(max(v, D.lbw(9)), D.ubw(9));
if abs(v - D.cw(9)) <= 1e-12 * (D.ubw(9) - D.lbw(9)), xh(9) = D.cx(9); else, xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9)); end
v = z(10); tl = 1e-12 * (D.ubw(10) - D.lbw(10)); if v < D.lbw(10) - tl || v > D.ubw(10) + tl, error('rastrigin_c_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
v = min(max(v, D.lbw(10)), D.ubw(10));
if abs(v - D.cw(10)) <= 1e-12 * (D.ubw(10) - D.lbw(10)), xh(10) = D.cx(10); else, xh(10) = map_work_to_orig(v, D.lbw(10), D.ubw(10), D.lb(10), D.ub(10)); end
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
if ~iscell(xr) || numel(xr) ~= 10, error('rastrigin_c_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('rastrigin_c_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('rastrigin_c_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('rastrigin_c_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('rastrigin_c_10D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = rastrigin_orig(x)
%
% Purpose:
%
%    Function rastrigin is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -5.12 <= x <= 5.12 (Montaz et al. (2005))
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
dim = size(x,1);
f   = 10*dim + sum(x.^2-10.*cos(2*pi.*x),1);
end
%
% End of rastrigin.

