function PB = tenn_local_minima_k_4D()
%TENN_LOCAL_MINIMA_K_4D  DFO Benchmark Suite v2.1.0 core: tenn_local_minima, dimension 4, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.177908578839979e-31 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K x3:K x4:K
%
lb_orig = [-10 -10 -10 -10];
ub_orig = [6.9961669921874998 6.4109028095310672 6.4787327999999995 6.9800083298625566];
lb_C_work = [-10 -10 -10 -10];
ub_C_work = [6.9961669921874998 6.4109028095310672 6.4787327999999995 6.9800083298625566];
T_break = [1.3307779947916671 0.94060187302071085 0.98582186666666693 1.320005553241705];   % breakpoints of T (NaN for continuous)
T_image = [1 1 1 1];   % their images
m = 25;
ProblemData = cell(1, 4);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = arrayfun(@(j) sprintf('k3_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'tenn_local_minima_k_4D';
PB.base_function = 'tenn_local_minima';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of tenn_local_minima (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'K','K','K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 4);
PB.proportion  = [0 0 4];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [25 25 25 25];
PB.n_finite    = 4;
PB.finite_part_size = 390625;
PB.M_K = 390625;
PB.domain_size = 390625;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'tenn_local_minima_k_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [NaN NaN NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([20 19 1 17 10 23 22 11 3 25 4 5 6 8 21 14 9 16 15 7 2 13 12 24 18])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([25 20 21 16 14 6 11 23 22 13 7 15 24 12 19 4 17 9 8 2 1 5 18 10 3])';
g = linspace(lb_orig(3), ub_orig(3), m); D.E{3} = g([24 3 16 4 12 17 22 7 6 14 9 11 15 2 13 10 8 23 21 19 1 18 25 5 20])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([10 15 19 5 18 7 25 22 8 13 21 2 23 6 3 17 24 11 14 16 1 9 4 20 12])';
D.fn = @(xh) tenn_local_minima_orig(xh);
PB.f_star      = 1.1779085788399792e-31;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [7 20 17 6];   % protocol starting point: Halton point 1058 (index representation)
PB.x0_alt_idx  = [14 21 11 7;1 22 22 5;25 20 15 6;12 21 3 7;18 22 20 4;6 21 8 6;22 21 25 7;9 22 13 8;15 20 1 5];   % 9 alternative starting points (Halton indices 1058 + 63 k)
PB.core_pidx   = 58;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('tenn_local_minima_k_4D: n = 4, pattern K, m = 25, f* = 1.177908579e-31 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('tenn_local_minima_k_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.E{3}(k);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
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
xr = cell(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 4, error('tenn_local_minima_k_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_k_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('tenn_local_minima_k_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('tenn_local_minima_k_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('tenn_local_minima_k_4D: x%d unknown label %s.', i, lab); end
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

