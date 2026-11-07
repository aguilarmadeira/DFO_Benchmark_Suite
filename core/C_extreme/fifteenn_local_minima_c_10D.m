function PB = fifteenn_local_minima_c_10D()
%FIFTEENN_LOCAL_MINIMA_C_10D  DFO Benchmark Suite v2 core: fifteenn_local_minima, dimension 10, pattern C, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10 10 10 10 10];
lb_C_work = [-10 -10 -10 -10 -10 -1000000000 -1000000000 -1000000000 -1000000000 -1000000000];
ub_C_work = [10 10 10 10 10 1000000000 1000000000 1000000000 1000000000 1000000000];
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
PB.name        = 'fifteenn_local_minima_c_10D';
PB.base_function = 'fifteenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of fifteenn_local_minima (F = f o T).';
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
D = struct('name', 'fifteenn_local_minima_c_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.fn = @(xh) fifteenn_local_minima_orig(xh);
PB.f_star      = 1.3497838043956716e-32;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [7.48046875 5.0983081847279372 8.5024000000000015 3.8359017076218258 3.5086401202103694 -225307237.14155662 938937512.72135162 346843563.20163298 -311744883.70181644 -709787199.14715648];   % protocol starting point: Halton point 1019 (index representation)
PB.x0_alt_idx  = [-2.802734375 3.6991312299954249 -0.61759999999999948 1.5118700541441061 -0.95416979714500272 -473827947.2007283 378383879.5033586 995043009.1850121 -908934001.8081696 -360203370.371889;2.353515625 4.2752629172382246 -6.6656000000000004 2.3865056226572268 -7.2201352366641629 -722348657.25989985 -182169753.71463478 -456188948.82635963 576888304.43001556 -10619541.596621513;-7.724609375 5.5098308184727927 3.4144000000000023 3.2611411911703456 8.1517655897821193 864360491.57942677 -859963362.50763273 192010497.15702009 66491329.004684687 338964287.17864633;8.681640625 3.8911751257430272 -2.7872000000000003 4.1357767596834645 1.8858001502629591 615839781.52025509 572562589.04946065 840209943.14039922 -443905646.42064595 688548115.95391369;-1.513671875 4.4673068129858269 7.2927999999999997 2.2115785089546041 -4.5440885185438162 367319071.46108341 12008955.831467509 -506050444.67123497 -954302621.84597683 -959489933.9866333;3.642578125 5.0434385002286231 1.2448000000000015 2.6863806747188654 -9.0068984358991866 -34137460.172963142 -548544677.38652563 142149001.31214476 531519684.39220834 -609906105.21136582;-6.435546875 3.3516232281664387 -7.8751999999999995 3.5610162432319843 4.5618468683833076 -282658170.2321347 766741298.59556293 685376877.09578633 21122708.96687746 -260322276.43609822;6.220703125 4.5861911294010067 5.2768000000000015 1.6368179925031221 0.099036951027937192 -684114701.86618114 206187665.37756991 -660883510.71584773 -489274266.45845318 89261552.339169502;-4.033203125 5.1623228166438047 -3.8431999999999995 2.511453561016241 -6.166928488491223 -932635411.92535269 -354365967.84042335 -12684064.732468367 996548039.77973199 438845381.1144371];   % 9 alternative starting points (Halton indices 1019 + 63 k)
PB.core_pidx   = 19;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('fifteenn_local_minima_c_10D: n = 10, pattern C, m = 0, f* = 1.349783804e-32 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('fifteenn_local_minima_c_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('fifteenn_local_minima_c_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('fifteenn_local_minima_c_10D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('fifteenn_local_minima_c_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('fifteenn_local_minima_c_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('fifteenn_local_minima_c_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
v = z(6); if v < D.lbw(6) || v > D.ubw(6), error('fifteenn_local_minima_c_10D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('fifteenn_local_minima_c_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
v = z(8); if v < D.lbw(8) || v > D.ubw(8), error('fifteenn_local_minima_c_10D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('fifteenn_local_minima_c_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
v = z(10); if v < D.lbw(10) || v > D.ubw(10), error('fifteenn_local_minima_c_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
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
if ~iscell(xr) || numel(xr) ~= 10, error('fifteenn_local_minima_c_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('fifteenn_local_minima_c_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('fifteenn_local_minima_c_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('fifteenn_local_minima_c_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('fifteenn_local_minima_c_10D: x%d unknown label %s.', i, lab); end
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

