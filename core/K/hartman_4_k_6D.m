function PB = hartman_4_k_6D()
%HARTMAN_4_K_6D  DFO Benchmark Suite v2.1.0 core: hartman_4, dimension 6, pattern K, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 17 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -3.322368011391339 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:K x2:K x3:K x4:K x5:K x6:K
%
lb_orig = [0 0.14795966062084032 0.14470336 0 0.18109077248821803 0.16527817653443508];
ub_orig = [0.84844116210937504 1 1 0.816347355268638 1 1];
lb_C_work = [0 0.14795966062084032 0.14470336 0 0.18109077248821803 0.16527817653443508];
ub_C_work = [0.84844116210937504 1 1 0.816347355268638 1 1];
T_break = [0.21211029052734376 0.14795966062084032 0.46543960000000001 0.25510854852144937 0.33463625264667718 0.63480920223381532];   % breakpoints of T (NaN for continuous)
T_image = [0.20169000000000001 0.15001100000000001 0.47687400000000002 0.27533200000000002 0.31165199999999998 0.6573];   % their images
m = 17;
ProblemData = cell(1, 6);
ProblemData{1} = arrayfun(@(j) sprintf('k1_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = arrayfun(@(j) sprintf('k3_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = arrayfun(@(j) sprintf('k5_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'hartman_4_k_6D';
PB.base_function = 'hartman_4';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: K adaptation of hartman_4 (F = f o T).';
PB.nvar        = 6;
PB.nobj        = 1;
PB.var_types   = {'K','K','K','K','K','K'};
PB.type_counts = struct('C', 0, 'D', 0, 'K', 6);
PB.proportion  = [0 0 6];
PB.pattern     = 'K';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [17 17 17 17 17 17];
PB.n_finite    = 6;
PB.finite_part_size = 24137569;
PB.M_K = 24137569;
PB.domain_size = 24137569;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'hartman_4_k_6D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 6);
D.cw = [NaN NaN NaN NaN NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN NaN NaN NaN NaN];
g = linspace(lb_orig(1), ub_orig(1), m); D.E{1} = g([17 3 9 8 15 14 16 12 5 6 2 7 11 4 10 1 13])';
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([7 13 4 8 11 6 17 2 10 1 9 5 15 16 14 12 3])';
g = linspace(lb_orig(3), ub_orig(3), m); D.E{3} = g([15 8 4 17 10 2 6 16 11 7 14 3 12 5 1 9 13])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([7 2 5 15 10 16 11 17 9 3 14 1 12 13 6 4 8])';
g = linspace(lb_orig(5), ub_orig(5), m); D.E{5} = g([13 2 15 1 6 10 4 3 5 9 14 17 12 8 7 16 11])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([17 15 13 14 7 1 8 16 9 4 5 10 12 2 11 6 3])';
D.fn = @(xh) hartman_4_orig(xh);
PB.f_star      = -3.322368011391339;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [5 1 4 12 6 16];   % protocol starting point: Halton point 1026 (index representation)
PB.x0_alt_idx  = [9 2 16 11 1 14;1 2 8 11 13 12;17 1 3 12 8 9;9 2 11 10 5 7;13 2 7 11 16 5;5 1 14 11 12 1;15 2 10 12 7 15;7 3 2 11 3 12;11 1 13 11 14 10];   % 9 alternative starting points (Halton indices 1026 + 63 k)
PB.core_pidx   = 26;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('hartman_4_k_6D: n = 6, pattern K, m = 17, f* = -3.322368011 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 6 || ~all(isfinite(z)), error('hartman_4_k_6D: z must be a finite numeric vector with 6 entries.'); end
xh = zeros(6, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.E{1}(k);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.E{3}(k);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
k = check_idx(z(5), D.m, D.name, 5); g = D.E{5}(k);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
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
xr = cell(1, 6);
for i = 1:6
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 6, error('hartman_4_k_6D: xr must be a cell array with 6 entries.'); end
z = zeros(1, 6);
for i = 1:6
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('hartman_4_k_6D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('hartman_4_k_6D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('hartman_4_k_6D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('hartman_4_k_6D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = hartman_4_orig(x)
%
% Purpose:
%
%    Function hartman_4 is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value
%    Global minima (one)
%    Local minima (4 local minima at x = A(:,n+2:2n+1))
%    Search domain: 0 <= x <= 1 (Brachetti et al. (1997))
%                               (Huyer and Neumaier (1999))
%                               (Huyer and Neumaier (2008))
%    Cases considered: n = 3 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      n = 6 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
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
if n == 3
    A = [3 10 30 1 0.3689 0.1170 0.2673
        0.1 10 35 1.2 0.4699 0.4387 0.7470
        3 10 30 3 0.1091 0.8732 0.5547
        0.1 10 35 3.2 0.03815 0.5743 0.8828];
else
    if n== 6
        A = [10 3 17 3.5 1.7 8 1 0.1312 0.1696 0.5569 0.0124 0.8283 0.5886
            0.05 10 17 0.1 8 14 1.2 0.2329 0.4135 0.8307 0.3736 0.1004 0.9991
            3 3.5 1.7 10 17 8 3 0.2348 0.1451 0.3522 0.2883 0.3047 0.6650
            17 8 0.05 10 0.1 14 3.2 0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
    end
end
f = -sum(A(:,n+1).*exp(-sum(A(:,1:n).*(repmat(x',4,1)-A(:,n+2:2*n+1)).^2,2)),1);
end
%
% End of hartman_4.

