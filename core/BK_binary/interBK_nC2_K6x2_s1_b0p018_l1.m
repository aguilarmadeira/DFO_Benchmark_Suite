function PB = interBK_nC2_K6x2_s1_b0p018_l1()
%INTERBK_NC2_K6X2_S1_B0P018_L1  DFO Benchmark Suite v2 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 6 categorical variables with m = 2 categories each (M_K = 64 combinations)
% lambda = 1, sigma = 0.15, b = 0.018, q = lambda^2 sigma^2 / b = 1.25, j_min = Inf (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.21581058502197267 0.76647793452142365];
s_ref = [-0.14543958461545617;-0.03670595628618821];
l2s = {[1 0], [0 1], [1 0], [1 0], [1 0], [1 0]};
k_ref = [0 0 1 0 0 1];
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2'}, {'k4_1';'k4_2'}, {'k5_1';'k5_2'}, {'k6_1';'k6_2'}, {'k7_1';'k7_2'}, {'k8_1';'k8_2'}};
D = struct('name', 'interBK_nC2_K6x2_s1_b0p018_l1', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 1, 'b', 0.017999999999999999);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K6x2_s1_b0p018_l1';
PB.base_function = 'familyBK';
PB.family      = 'BK';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core, Family BK (controlled categorical dimension): 6 categorical variables x 2 categories.';
PB.nvar        = 8;
PB.nobj        = 1;
PB.var_types   = {'C','C','K','K','K','K','K','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 6);
PB.proportion  = [2 0 6];
PB.pattern     = 'BK';
PB.configuration = 'interBK_nC2_K6x2';
PB.nK = 6; PB.m = 2; PB.M_K = 64;
PB.lambda      = 1;
PB.sigma       = 0.14999999999999999;
PB.b           = 0.017999999999999999;
PB.q           = 1.25;
PB.j_min       = Inf;
PB.obstacle_ratio = 1.25;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 1;
PB.relabel     = 1;
PB.lb_orig     = zeros(8,1); PB.ub_orig = ones(8,1);
PB.lb_C_work   = PB.lb_orig; PB.ub_C_work = PB.ub_orig;
PB.ProblemData = ProblemData;
PB.cardinality = [Inf Inf [2 2 2 2 2 2]];
PB.cardinality = PB.cardinality(:)';
PB.n_finite    = 6;
PB.finite_part_size = 64;
PB.domain_size = Inf;
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by construction
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.921142578125 0.57262612406645341 1 1 2 1 1 2];   % protocol starting point: Halton point 3031 (index representation)
PB.x0_alt_idx  = [0.406982421875 0.6014327084285932 2 1 1 1 1 2;0.664794921875 0.66316110349032165 1 1 1 2 2 1;0.160888671875 0.5822283188538333 1 1 2 2 1 2;0.793701171875 0.6110349032159732 2 1 1 2 1 2;0.283935546875 0.6398414875781131 1 1 1 1 2 1;0.541748046875 0.55890870294162487 2 1 2 1 2 1;0.037841796875 0.6206370980033532 2 1 2 1 1 2;0.951904296875 0.6494436823654931 1 1 1 2 2 1;0.439208984375 0.57948483462886768 2 1 1 2 2 1];   % 9 alternative starting points (Halton indices 3031 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K6x2_s1_b0p018_l1: n = 8 (nC = 2, nK = 6 x 2), lambda = 1, b = 0.018, j_min = Inf, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 8 || ~all(isfinite(z)), error('interBK_nC2_K6x2_s1_b0p018_l1: z must be a finite numeric vector with 8 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K6x2_s1_b0p018_l1: continuous coordinates must lie in [0,1].'); end
kk = z(3:8);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 2), error('interBK_nC2_K6x2_s1_b0p018_l1: label indices must be integers in 1..2.'); end
h = 0;
for i = 1:6, h = h + (D.l2s{i}(kk(i)) == D.ks(i)); end
if h == 0, shift = D.lambda * D.s; else, shift = zeros(2,1); end
dd = u - D.xs - shift;
f = dd' * dd + D.b * (1 - h / 6);
end

function xr = idx2cell(z, PB)
xr = cell(1, 8);
for i = 1:2, xr{i} = z(i); end
for i = 3:8
    k = z(i); if k ~= floor(k) || k < 1 || k > 2, error('interBK_nC2_K6x2_s1_b0p018_l1: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 8, error('interBK_nC2_K6x2_s1_b0p018_l1: xr must be a cell array with 8 entries.'); end
z = zeros(1, 8);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K6x2_s1_b0p018_l1: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:8
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K6x2_s1_b0p018_l1: unknown label %s.', lab); end
    z(i) = j;
end
end
