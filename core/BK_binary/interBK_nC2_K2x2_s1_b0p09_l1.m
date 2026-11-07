function PB = interBK_nC2_K2x2_s1_b0p09_l1()
%INTERBK_NC2_K2X2_S1_B0P09_L1  DFO Benchmark Suite v2 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 2 categorical variables with m = 2 categories each (M_K = 4 combinations)
% lambda = 1, sigma = 0.15, b = 0.09, q = lambda^2 sigma^2 / b = 0.25, j_min = 1 (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.28789920806884767 0.29710504082297001];
s_ref = [-0.14502423761214869;-0.038314103212982002];
l2s = {[1 0], [1 0]};
k_ref = [1 1];
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2'}, {'k4_1';'k4_2'}};
D = struct('name', 'interBK_nC2_K2x2_s1_b0p09_l1', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 1, 'b', 0.089999999999999997);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K2x2_s1_b0p09_l1';
PB.base_function = 'familyBK';
PB.family      = 'BK';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core, Family BK (controlled categorical dimension): 2 categorical variables x 2 categories.';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','C','K','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'BK';
PB.configuration = 'interBK_nC2_K2x2';
PB.nK = 2; PB.m = 2; PB.M_K = 4;
PB.lambda      = 1;
PB.sigma       = 0.14999999999999999;
PB.b           = 0.089999999999999997;
PB.q           = 0.25;
PB.j_min       = 1;
PB.obstacle_ratio = 0.25;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 1;
PB.relabel     = 1;
PB.lb_orig     = zeros(4,1); PB.ub_orig = ones(4,1);
PB.lb_C_work   = PB.lb_orig; PB.ub_C_work = PB.ub_orig;
PB.ProblemData = ProblemData;
PB.cardinality = [Inf Inf [2 2]];
PB.cardinality = PB.cardinality(:)';
PB.n_finite    = 2;
PB.finite_part_size = 4;
PB.domain_size = Inf;
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by construction
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.764892578125 0.81953970431336676 1 1];   % protocol starting point: Halton point 3011 (index representation)
PB.x0_alt_idx  = [0.250732421875 0.8812680993750952 2 1;0.508544921875 0.81130925163846968 1 1;0.004638671875 0.82914189910074676 1 1;0.989013671875 0.85794848346288666 2 1;0.496826171875 0.78798963572626113 1 1;0.737060546875 0.84971803078798946 2 1;0.244873046875 0.86755067825026666 1 1;0.865966796875 0.79759183051364124 1 1;0.373779296875 0.82639841487578103 2 1];   % 9 alternative starting points (Halton indices 3011 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K2x2_s1_b0p09_l1: n = 4 (nC = 2, nK = 2 x 2), lambda = 1, b = 0.09, j_min = 1, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('interBK_nC2_K2x2_s1_b0p09_l1: z must be a finite numeric vector with 4 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K2x2_s1_b0p09_l1: continuous coordinates must lie in [0,1].'); end
kk = z(3:4);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 2), error('interBK_nC2_K2x2_s1_b0p09_l1: label indices must be integers in 1..2.'); end
h = 0;
for i = 1:2, h = h + (D.l2s{i}(kk(i)) == D.ks(i)); end
if h == 0, shift = D.lambda * D.s; else, shift = zeros(2,1); end
dd = u - D.xs - shift;
f = dd' * dd + D.b * (1 - h / 2);
end

function xr = idx2cell(z, PB)
xr = cell(1, 4);
for i = 1:2, xr{i} = z(i); end
for i = 3:4
    k = z(i); if k ~= floor(k) || k < 1 || k > 2, error('interBK_nC2_K2x2_s1_b0p09_l1: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 4, error('interBK_nC2_K2x2_s1_b0p09_l1: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K2x2_s1_b0p09_l1: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:4
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K2x2_s1_b0p09_l1: unknown label %s.', lab); end
    z(i) = j;
end
end
