function PB = fifteenn_local_minima_c_8D()
%FIFTEENN_LOCAL_MINIMA_C_8D  DFO Benchmark Suite v2 core: fifteenn_local_minima, dimension 8, pattern C, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10 10 10];
lb_C_work = [-10 -10 -10 -10 -10 -10 -10 -10];
ub_C_work = [10 10 10 10 10 10 10 10];
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
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of fifteenn_local_minima (F = f o T).';
PB.nvar        = 8;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 8, 'D', 0, 'K', 0);
PB.proportion  = [8 0 0];
PB.pattern     = 'C';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1 1 1 1 1 1];
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
D.fn = @(xh) fifteenn_local_minima_orig(xh);
PB.f_star      = 1.3497838043956716e-32;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-2.51953125 -1.5683584819387288 4.5024000000000015 0.97875885047896638 1.6904583020285493 -3.7915339098771046 8.212904538978222 2.4158040530689622];   % protocol starting point: Halton point 1018 (index representation)
PB.x0_alt_idx  = [2.197265625 -2.9675354366712394 -4.6175999999999995 -1.3452728029987515 -2.7723516153268228 -6.2767410104688208 2.6073682067982915 8.897798512902753;-7.646484375 -2.3914037494284397 8.534399999999998 -0.47063723448563266 -9.0383170548459795 -8.7619481110605371 -2.9981681253816408 -5.6145210672109638;8.525390625 -1.1568358481938716 -0.58559999999999945 0.40399833402748797 6.3335837716002992 7.1051433773327304 -9.7761042133116227 0.86747339262283063;-1.318359375 -2.7754915409236389 -6.7872000000000003 1.2786339025406068 0.06761833208114254 4.6199362767410097 4.5491553022593116 7.3494678524566268;3.486328125 -2.199359853680841 3.2928000000000033 -0.64556434818825537 -6.3622703367256328 2.1347291761492961 -1.0563810299206207 -6.1131360256597178;-6.357421875 -1.623228166438043 -2.7552000000000012 -0.17076218242399044 9.0096304897206494 -1.8798361401911698 -6.6619173621005494 0.36885843417407749;6.064453125 -3.3150434385002292 7.3248000000000033 0.70387338608912664 2.7436650502014892 -4.3650432407828852 6.4909423977203353 5.801137192010497;-3.779296875 -2.0804755372656611 1.2767999999999997 -1.2203248646397338 -1.7191448671538829 -8.3796085571233494 0.8854060655404048 -7.6614666861058467;0.966796875 -1.5043438500228632 -7.8431999999999995 -0.34568929612661492 -7.9851103066730413 9.0168411470186633 -4.7201302666395284 -1.1794722262720523];   % 9 alternative starting points (Halton indices 1018 + 63 k)
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
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('fifteenn_local_minima_c_8D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('fifteenn_local_minima_c_8D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('fifteenn_local_minima_c_8D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('fifteenn_local_minima_c_8D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('fifteenn_local_minima_c_8D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
v = z(6); if v < D.lbw(6) || v > D.ubw(6), error('fifteenn_local_minima_c_8D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('fifteenn_local_minima_c_8D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
v = z(8); if v < D.lbw(8) || v > D.ubw(8), error('fifteenn_local_minima_c_8D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8));
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
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('fifteenn_local_minima_c_8D: x%d must be a grid value.', i); end
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

