function PB = transistor_ck_9D()
%TRANSISTOR_CK_9D  DFO Benchmark Suite v2 core: transistor, dimension 9, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 25 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C x6:K x7:C x8:K x9:C
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10 10 10 10];
lb_C_work = [-10 -10 -10 -10 -10 -10 -10 -10 -10];
ub_C_work = [10 10 10 10 10 10 10 10 10];
T_break = [NaN 0.83333333333333337 NaN 1.6666666666666667 NaN 8.3333333333333339 NaN 0.83333333333333337 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0.44998747198158512 NaN 2.000068541623683 NaN 7.9996926842185427 NaN 0.99998772345678977 NaN];   % their images
m = 25;
ProblemData = cell(1, 9);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = arrayfun(@(j) sprintf('k6_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = arrayfun(@(j) sprintf('k8_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
PB.name        = 'transistor_ck_9D';
PB.base_function = 'transistor';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of transistor (F = f o T).';
PB.nvar        = 9;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C','K','C','K','C'};
PB.type_counts = struct('C', 5, 'D', 0, 'K', 4);
PB.proportion  = [5 0 4];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1 1 1];
PB.cardinality = [Inf 25 Inf 25 Inf 25 Inf 25 Inf];
PB.n_finite    = 4;
PB.finite_part_size = 390625;
PB.M_K = 390625;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'transistor_ck_9D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 9);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([2 5 22 21 12 11 25 8 17 13 10 24 3 20 1 9 15 16 4 19 23 14 7 18 6])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([2 23 17 9 25 5 22 12 15 21 4 3 19 20 11 6 14 1 10 13 24 18 16 7 8])';
g = linspace(lb_orig(6), ub_orig(6), m); D.E{6} = g([17 11 5 16 13 24 3 6 8 1 12 9 18 19 22 23 14 4 15 21 7 20 25 10 2])';
g = linspace(lb_orig(8), ub_orig(8), m); D.E{8} = g([18 9 4 6 21 17 10 16 24 5 25 15 23 1 3 14 13 7 20 22 8 2 11 19 12])';
D.fn = @(xh) transistor_orig(xh);
PB.f_star      = 1.6004932665473596e-24;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-2.177734375 2 0.022400000000001086 20 2.3516153268219391 18 0.18522287807856763 24 -6.5184515492726227];   % protocol starting point: Halton point 1062 (index representation)
PB.x0_alt_idx  = [2.978515625 4 -9.8656000000000006 21 -3.9143501126972202 15 -5.4203134541013638 6 8.3397715131092269;-7.099609375 2 4.0544000000000011 19 -8.3771600300525915 10 7.7325463057195201 14 3.2358017588559207;5.556640625 2 -5.9871999999999996 20 5.1915852742299045 7 2.1270099735395913 22 -1.8681679953973873;-4.638671875 3 7.9328000000000003 21 0.72877535687453232 4 -3.4785263586403428 6 -6.9721377496506944;0.517578125 1 -1.9551999999999996 18 -5.7011133119322448 23 -9.0840626908202733 13 7.8860853127311579;-9.560546875 3 -7.2351999999999999 19 9.6707875145140356 20 5.2411968247506593 21 2.7821155584778481;9.658203125 3 2.0768000000000004 20 3.4048220749948772 17 -1.5367392631793209 4 -2.3218541957754582;-0.595703125 1 -3.2032000000000007 21 -1.0579878423604949 13 -7.1422755953592514 12 -7.4258239500287662;4.560546875 2 6.1088000000000022 19 -7.3239532818796524 10 7.182983920211683 19 6.5644776855428617];   % 9 alternative starting points (Halton indices 1062 + 63 k)
PB.core_pidx   = 62;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('transistor_ck_9D: n = 9, pattern CK, m = 25, f* = 1.600493267e-24 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 9 || ~all(isfinite(z)), error('transistor_ck_9D: z must be a finite numeric vector with 9 entries.'); end
xh = zeros(9, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('transistor_ck_9D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('transistor_ck_9D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('transistor_ck_9D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
k = check_idx(z(6), D.m, D.name, 6); g = D.E{6}(k);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('transistor_ck_9D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
k = check_idx(z(8), D.m, D.name, 8); g = D.E{8}(k);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('transistor_ck_9D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
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
xr = cell(1, 9);
for i = 1:9
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 9, error('transistor_ck_9D: xr must be a cell array with 9 entries.'); end
z = zeros(1, 9);
for i = 1:9
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('transistor_ck_9D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('transistor_ck_9D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('transistor_ck_9D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('transistor_ck_9D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = transistor_orig(x)
%
% Purpose:
%
%    Function transistor is the function described in
%    Montaz et al. (2005).
%
%    dim = 9
%    Minimum global value = 0
%    Global minima near [0.9 0.45 1 2 8 8 5 1 2]'
%    Local minima
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
g = [0.485 0.752 0.869 0.982
    0.369 1.254 0.703 1.455
    5.2095 10.0677 22.9274 20.2153
    23.3037 101.779 111.461 191.267
    28.5132 111.8467 134.3884 211.4823]';
alfa = (1-x(1)*x(2))*x(3).*(exp(x(5).*(g(:,1)-g(:,3).*x(7)*0.001-...
       g(:,5).*x(8)*0.001))-1)-g(:,5)+g(:,4).*x(2);
beta = (1-x(1)*x(2))*x(4).*(exp(x(6).*(g(:,1)-g(:,2)-g(:,3).*x(7)*0.001+...
       g(:,4).*x(9)*0.001))-1)-g(:,5).*x(1)+g(:,4);
f = (x(1)*x(3)-x(2)*x(4))^2 + sum(alfa.^2+beta.^2);
end
%
% End of transistor.

