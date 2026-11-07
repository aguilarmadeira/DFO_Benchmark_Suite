function PB = fifteenn_local_minima_cdk_8D()
%FIFTEENN_LOCAL_MINIMA_CDK_8D  DFO Benchmark Suite v2 core: fifteenn_local_minima, dimension 8, pattern CDK, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.349783804395672e-32 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C x5:D x6:K x7:C x8:K
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10 10 10];
lb_C_work = [-10 -10 -10 -1000000000 -10 -10 -1000000000 -10];
ub_C_work = [10 10 10 1000000000 10 10 1000000000 10];
T_break = [NaN 1.25 1.25 NaN 1.25 1.25 NaN 1.25];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 1 NaN 1 1 NaN 1];   % their images
m = 17;
ProblemData = cell(1, 8);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = {linspace(lb_orig(5), ub_orig(5), m)'};
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'fifteenn_local_minima_cdk_8D';
PB.base_function = 'fifteenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CDK adaptation of fifteenn_local_minima (F = f o T).';
PB.nvar        = 8;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C','D','K','C','K'};
PB.type_counts = struct('C', 3, 'D', 2, 'K', 3);
PB.proportion  = [3 2 3];
PB.pattern     = 'CDK';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 100000000 100000000];
PB.cardinality = [Inf 17 17 Inf 17 17 Inf 17];
PB.n_finite    = 5;
PB.finite_part_size = 1419857;
PB.M_K = 4913;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'fifteenn_local_minima_cdk_8D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 8);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([7 12 14 15 17 6 9 4 10 3 16 8 5 13 11 2 1])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([5 13 14 11 7 12 10 16 2 6 9 1 17 3 8 15 4])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([11 3 12 10 7 13 4 15 16 8 1 9 5 17 14 2 6])';
D.fn = @(xh) fifteenn_local_minima_orig(xh);
PB.f_star      = 1.3497838043956716e-32;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-2.51953125 8 13 97875885.047896624 10 6 821290453.89782214 11];   % protocol starting point: Halton point 1018 (index representation)
PB.x0_alt_idx  = [2.197265625 7 5 -134527280.29987514 7 4 260736820.67982912 16;-7.646484375 7 16 -47063723.448563218 2 2 -299816812.53816402 5;8.525390625 8 9 40399833.402748942 14 15 -977610421.33116221 10;-1.318359375 7 4 127863390.25406075 9 13 454915530.22593117 15;3.486328125 7 12 -64556434.818825603 4 11 -105638102.99206197 4;-6.357421875 8 7 -17076218.242398977 16 7 -666191736.21005487 9;6.064453125 6 15 70387338.608912706 11 6 649094239.77203345 14;-3.779296875 7 10 -122032486.4639734 8 2 88540606.554040432 3;0.966796875 8 3 -34568929.612661481 3 16 -472013026.66395283 8];   % 9 alternative starting points (Halton indices 1018 + 63 k)
PB.core_pidx   = 18;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('fifteenn_local_minima_cdk_8D: n = 8, pattern CDK, m = 17, f* = 1.349783804e-32 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 8 || ~all(isfinite(z)), error('fifteenn_local_minima_cdk_8D: z must be a finite numeric vector with 8 entries.'); end
xh = zeros(8, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('fifteenn_local_minima_cdk_8D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('fifteenn_local_minima_cdk_8D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
k = check_idx(z(5), D.m, D.name, 5); g = D.lb(5) + (k - 1) * (D.ub(5) - D.lb(5)) / (D.m - 1);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('fifteenn_local_minima_cdk_8D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
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
if ~iscell(xr) || numel(xr) ~= 8, error('fifteenn_local_minima_cdk_8D: xr must be a cell array with 8 entries.'); end
z = zeros(1, 8);
for i = 1:8
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('fifteenn_local_minima_cdk_8D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('fifteenn_local_minima_cdk_8D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('fifteenn_local_minima_cdk_8D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('fifteenn_local_minima_cdk_8D: x%d unknown label %s.', i, lab); end
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

