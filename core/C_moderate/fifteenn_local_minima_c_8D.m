function PB = fifteenn_local_minima_c_8D()
%FIFTEENN_LOCAL_MINIMA_C_8D  DFO Benchmark Suite v2.1.0 core: fifteenn_local_minima, dimension 8, pattern C, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [6.9250732421874996 6.4877203678301072 6.3891327999999996 6.5555185339441904 7.3963663684174579 6.1642589545183988 6.225670190730475 6.8531748528633143];
lb_C_work = [-10 -1000 -0.10000000000000001 -10 -1000 -0.10000000000000001 -10 -1000];
ub_C_work = [6.9250732421875005 648.77203678301066 0.063891327999999997 6.5555185339441913 739.63663684174571 0.061642589545183986 6.2256701907304759 685.31748528633148];
T_break = [NaN NaN NaN NaN NaN NaN NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN NaN NaN NaN NaN NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 8);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = [lb_C_work(6); ub_C_work(6)];
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = [lb_C_work(8); ub_C_work(8)];
PB.name        = 'fifteenn_local_minima_c_8D';
PB.base_function = 'fifteenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: C adaptation of fifteenn_local_minima (F = f o T).';
PB.nvar        = 8;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 8, 'D', 0, 'K', 0);
PB.proportion  = [8 0 0];
PB.pattern     = 'C';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01 1 100 0.01 1 100];
PB.cardinality = [Inf Inf Inf Inf Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'fifteenn_local_minima_c_8D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 8);
D.cw = [1 100 0.009999999999999995 1 100 0.009999999999999995 1 100];   % anchors of the continuous maps (NaN for finite)
D.cx = [1 1 1 1 1 1 1 1];
D.fn = @(xh) fifteenn_local_minima_orig(xh);
PB.f_star      = 1.3497838043956716e-32;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-3.6696259260177611 -304.90726204209614 0.01884087975936001 -0.91204771855958278 16.857478183980561 -0.049822373204453589 4.7758291132359343 46.228583226301225];   % protocol starting point: Halton point 1018 (index representation)
PB.x0_alt_idx  = [0.32198070287704539 -420.2534539158047 -0.055893565808640001 -2.8358251741907488 -371.32590359062385 -0.069908138769239764 0.22814992483050744 592.43951336065197;-8.0083287835121162 -372.7579631442776 0.051881371484159999 -2.1118229059424607 -916.35105577920979 -0.089993904334025931 -4.3195292635749203 -630.4537836607899;5.6771796584129355 -270.98191149100535 -0.022853074083840005 -1.3878206376941709 420.7250369999831 0.038245983502685821 -9.8183570404050116 -84.242853526438694;-2.6531298279762261 -404.42162365862896 -0.073672497070080001 -0.66381836944588457 -124.30011518860272 0.018160217937899611 1.8034897744088578 461.96807660791251;1.412854564189912 -356.92613288710186 0.0089287322419200205 -2.2566233595921181 -683.58361014209618 -0.0019255476268865285 -2.7441894139965699 -672.47000905573998;-6.9174549221992496 -309.43064211557476 -0.040632005345280012 -1.8635935568287607 653.49248263709683 -0.034371784308464229 -7.2918686024019967 -126.25907892138878;3.5946022868156451 -448.90152771450357 0.04196922396672001 -1.1395912885804744 108.46733044851067 -0.054457549873250397 3.37882962398721 331.49663985517282;-4.7357071995735165 -347.12547606123132 -0.0075915136204800088 -2.7323962787267062 -279.71605132609375 -0.086903786554828069 -1.1688495644182186 -802.94144580847967;-0.71930798292159892 -299.62998528970434 -0.08232595918848 -2.008394010478419 -824.74120351467968 0.053696572398675174 -5.7165287528236464 -256.73051567412858];   % 9 alternative starting points (Halton indices 1018 + 63 k)
PB.core_pidx   = 18;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('fifteenn_local_minima_c_8D: n = 8, pattern C, m = 0, f* = 1.349783804e-32 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 8 || ~all(isfinite(z)), error('fifteenn_local_minima_c_8D: z must be a finite numeric vector with 8 entries.'); end
xh = zeros(8, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('fifteenn_local_minima_c_8D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
v = z(2); tl = 1e-12 * (D.ubw(2) - D.lbw(2)); if v < D.lbw(2) - tl || v > D.ubw(2) + tl, error('fifteenn_local_minima_c_8D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
v = min(max(v, D.lbw(2)), D.ubw(2));
if abs(v - D.cw(2)) <= 1e-12 * (D.ubw(2) - D.lbw(2)), xh(2) = D.cx(2); else, xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2)); end
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('fifteenn_local_minima_c_8D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('fifteenn_local_minima_c_8D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('fifteenn_local_minima_c_8D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
v = z(6); tl = 1e-12 * (D.ubw(6) - D.lbw(6)); if v < D.lbw(6) - tl || v > D.ubw(6) + tl, error('fifteenn_local_minima_c_8D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
v = min(max(v, D.lbw(6)), D.ubw(6));
if abs(v - D.cw(6)) <= 1e-12 * (D.ubw(6) - D.lbw(6)), xh(6) = D.cx(6); else, xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6)); end
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('fifteenn_local_minima_c_8D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
v = z(8); tl = 1e-12 * (D.ubw(8) - D.lbw(8)); if v < D.lbw(8) - tl || v > D.ubw(8) + tl, error('fifteenn_local_minima_c_8D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
v = min(max(v, D.lbw(8)), D.ubw(8));
if abs(v - D.cw(8)) <= 1e-12 * (D.ubw(8) - D.lbw(8)), xh(8) = D.cx(8); else, xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8)); end
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
xr = cell(1, 8);
for i = 1:8
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 8, error('fifteenn_local_minima_c_8D: xr must be a cell array with 8 entries.'); end
z = zeros(1, 8);
for i = 1:8
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('fifteenn_local_minima_c_8D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('fifteenn_local_minima_c_8D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('fifteenn_local_minima_c_8D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('fifteenn_local_minima_c_8D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = fifteenn_local_minima_orig(x)
%
% Purpose:
%
%    Function fifteenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (15^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
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
f   = 0.1*sum(((x(1:dim-1)-1).^2).*(1+10*sin(3*pi.*x(2:dim)).^2),1);
f   = f + 0.1*(sin(3*pi*x(1))^2+(x(dim)-1)^2*(1+sin(2*pi*x(dim))^2));
end
%
% End of fifteenn_local_minima.

