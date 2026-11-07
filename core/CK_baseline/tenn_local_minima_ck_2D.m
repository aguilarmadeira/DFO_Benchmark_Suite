function PB = tenn_local_minima_ck_2D()
%TENN_LOCAL_MINIMA_CK_2D  DFO Benchmark Suite v2 core: tenn_local_minima, dimension 2, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 65 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 2.355817157679958e-31 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K
%
lb_orig = [-10 -10];
ub_orig = [10 10];
lb_C_work = [-10 -10];
ub_C_work = [10 10];
T_break = [NaN 0.9375];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1];   % their images
m = 65;
ProblemData = cell(1, 2);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'tenn_local_minima_ck_2D';
PB.base_function = 'tenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of tenn_local_minima (F = f o T).';
PB.nvar        = 2;
PB.nobj        = 1;
PB.var_types   = {'C','K'};
PB.type_counts = struct('C', 1, 'D', 0, 'K', 1);
PB.proportion  = [1 0 1];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = 1;
PB.cardinality = [Inf 65];
PB.n_finite    = 1;
PB.finite_part_size = 65;
PB.M_K = 65;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'tenn_local_minima_ck_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([60 17 55 7 54 4 65 44 36 16 46 34 11 35 49 2 9 47 58 33 12 63 29 8 38 31 50 52 43 27 51 37 1 59 3 14 30 42 18 53 32 24 6 57 25 28 23 56 41 5 26 21 19 10 62 48 45 20 40 39 64 22 15 13 61])';
D.fn = @(xh) tenn_local_minima_orig(xh);
PB.f_star      = 2.3558171576799584e-31;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.322265625 30];   % protocol starting point: Halton point 1057 (index representation)
PB.x0_alt_idx  = [-9.521484375 34;9.462890625 36;-0.380859375 30;4.423828125 32;-5.419921875 36;7.001953125 32;-2.841796875 33;1.904296875 35;-7.939453125 30];   % 9 alternative starting points (Halton indices 1057 + 63 k)
PB.core_pidx   = 57;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('tenn_local_minima_ck_2D: n = 2, pattern CK, m = 65, f* = 2.355817158e-31 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('tenn_local_minima_ck_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('tenn_local_minima_ck_2D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
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
xr = cell(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 2, error('tenn_local_minima_ck_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_ck_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('tenn_local_minima_ck_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('tenn_local_minima_ck_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('tenn_local_minima_ck_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = tenn_local_minima_orig(x)
%
% Purpose:
%
%    Function tenn_local_minima is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (10^n different points)
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
%                      n = 6 Brachetti et al. (1997)
%                      n = 8 Brachetti et al. (1997)
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
f   = pi/dim*sum(((x(1:dim-1)-1).^2).*(1+10*sin(pi.*x(2:dim)).^2),1);
f   = f + pi/dim*(10*sin(pi*x(1))^2+(x(dim)-1)^2);
end
%
% End of tenn_local_minima.

