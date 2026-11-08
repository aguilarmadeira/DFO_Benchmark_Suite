function PB = hartman_4_dk_3D()
%HARTMAN_4_DK_3D  DFO Benchmark Suite v2.1.0 core: hartman_4, dimension 3, pattern DK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -3.862782147819745 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:K x3:D
%
lb_orig = [0 0 0];
ub_orig = [0.83094116210937496 0.81315145049027082 0.86929663999999995];
lb_C_work = [0 0 0];
ub_C_work = [0.83094116210937496 0.81315145049027082 0.86929663999999995];
T_break = [0.12464117431640624 0.54887722908093273 0.84756422399999987];   % breakpoints of T (NaN for continuous)
T_image = [0.11461399999999999 0.55564899999999995 0.85254700000000005];   % their images
m = 41;
ProblemData = cell(1, 3);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
PB.name        = 'hartman_4_dk_3D';
PB.base_function = 'hartman_4';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: DK adaptation of hartman_4 (F = f o T).';
PB.nvar        = 3;
PB.nobj        = 1;
PB.var_types   = {'D','K','D'};
PB.type_counts = struct('C', 0, 'D', 2, 'K', 1);
PB.proportion  = [0 2 1];
PB.pattern     = 'DK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [41 41 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 41;
PB.domain_size = 68921;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'hartman_4_dk_3D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 3);
D.cw = [NaN NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([10 21 40 16 37 24 7 18 39 17 36 11 15 30 5 12 29 32 33 25 1 20 35 2 28 14 13 6 19 9 8 27 26 34 4 31 3 41 23 22 38])';
D.fn = @(xh) hartman_4_orig(xh);
PB.f_star      = -3.862782147819745;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [21 40 2];   % protocol starting point: Halton point 1025 (index representation)
PB.x0_alt_idx  = [1 37 29;41 38 10;21 41 37;31 38 17;10 39 7;36 40 25;16 37 15;26 39 33;5 40 23];   % 9 alternative starting points (Halton indices 1025 + 63 k)
PB.core_pidx   = 25;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('hartman_4_dk_3D: n = 3, pattern DK, m = 41, f* = -3.862782148 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('hartman_4_dk_3D: z must be a finite numeric vector with 3 entries.'); end
xh = zeros(3, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.lb(1) + (k - 1) * (D.ub(1) - D.lb(1)) / (D.m - 1);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
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
xr = cell(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 3, error('hartman_4_dk_3D: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('hartman_4_dk_3D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('hartman_4_dk_3D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('hartman_4_dk_3D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('hartman_4_dk_3D: x%d unknown label %s.', i, lab); end
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

