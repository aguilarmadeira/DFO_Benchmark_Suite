function PB = hartman_4_cd_6D()
%HARTMAN_4_CD_6D  DFO Benchmark Suite v2 core: hartman_4, dimension 6, pattern CD, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -3.322368011391339 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:D x3:C x4:D x5:C x6:D
%
lb_orig = [0 0 0 0 0 0];
ub_orig = [1 1 1 1 1 1];
lb_C_work = [0 0 0 0 0 0];
ub_C_work = [1 1 1 1 1 1];
T_break = [NaN 0.14999999999999999 NaN 0.27500000000000002 NaN 0.65000000000000002];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0.15001100000000001 NaN 0.27533200000000002 NaN 0.6573];   % their images
m = 41;
ProblemData = cell(1, 6);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = {linspace(lb_orig(2), ub_orig(2), m)'};
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = {linspace(lb_orig(4), ub_orig(4), m)'};
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = {linspace(lb_orig(6), ub_orig(6), m)'};
PB.name        = 'hartman_4_cd_6D';
PB.base_function = 'hartman_4';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CD adaptation of hartman_4 (F = f o T).';
PB.nvar        = 6;
PB.nobj        = 1;
PB.var_types   = {'C','D','C','D','C','D'};
PB.type_counts = struct('C', 3, 'D', 3, 'K', 0);
PB.proportion  = [3 3 0];
PB.pattern     = 'CD';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1];
PB.cardinality = [Inf 41 Inf 41 Inf 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'hartman_4_cd_6D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 6);
D.fn = @(xh) hartman_4_orig(xh);
PB.f_star      = -3.322368011391339;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.25048828125 2 0.21312 29 0.32006010518407207 38];   % protocol starting point: Halton point 1026 (index representation)
PB.x0_alt_idx  = [0.50830078125 3 0.90912000000000004 25 0.0067618332081142004 33;0.00439453125 4 0.41472000000000003 27 0.77535687453042823 28;0.98876953125 1 0.15071999999999997 28 0.46205860255447034 20;0.49658203125 4 0.60864000000000007 24 0.23891810668670174 15;0.73681640625 5 0.34464 26 0.90915921043644565 10;0.24462890625 2 0.81024000000000007 27 0.68601871456867691 2;0.86572265625 3 0.54624000000000006 29 0.37272044259271903 37;0.37353515625 5 0.051840000000000004 25 0.14957994672495048 29;0.61083984375 2 0.74784000000000006 27 0.82801721193907518 24];   % 9 alternative starting points (Halton indices 1026 + 63 k)
PB.core_pidx   = 26;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('hartman_4_cd_6D: n = 6, pattern CD, m = 41, f* = -3.322368011 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 6 || ~all(isfinite(z)), error('hartman_4_cd_6D: z must be a finite numeric vector with 6 entries.'); end
xh = zeros(6, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('hartman_4_cd_6D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('hartman_4_cd_6D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.lb(4) + (k - 1) * (D.ub(4) - D.lb(4)) / (D.m - 1);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('hartman_4_cd_6D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.lb(6) + (k - 1) * (D.ub(6) - D.lb(6)) / (D.m - 1);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
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
if ~iscell(xr) || numel(xr) ~= 6, error('hartman_4_cd_6D: xr must be a cell array with 6 entries.'); end
z = zeros(1, 6);
for i = 1:6
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('hartman_4_cd_6D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('hartman_4_cd_6D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('hartman_4_cd_6D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('hartman_4_cd_6D: x%d unknown label %s.', i, lab); end
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

