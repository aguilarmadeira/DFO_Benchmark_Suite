function PB = interBK_nC2_K4x2_s3_b0p06_l0()
%INTERBK_NC2_K4X2_S3_B0P06_L0  DFO Benchmark Suite v2 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 4 categorical variables with m = 2 categories each (M_K = 16 combinations)
% lambda = 0, sigma = 0.15, b = 0.06, q = lambda^2 sigma^2 / b = 0, j_min = 0 (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.29114017486572263 0.32994556310107803];
s_ref = [-0.056297606711978164;-0.13903445428563177];
l2s = {[1 0], [1 0], [1 0], [1 0]};
k_ref = [1 0 1 1];
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2'}, {'k4_1';'k4_2'}, {'k5_1';'k5_2'}, {'k6_1';'k6_2'}};
D = struct('name', 'interBK_nC2_K4x2_s3_b0p06_l0', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 0, 'b', 0.059999999999999998);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K4x2_s3_b0p06_l0';
PB.base_function = 'familyBK';
PB.family      = 'BK';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core, Family BK (controlled categorical dimension): 4 categorical variables x 2 categories.';
PB.nvar        = 6;
PB.nobj        = 1;
PB.var_types   = {'C','C','K','K','K','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 4);
PB.proportion  = [2 0 4];
PB.pattern     = 'BK';
PB.configuration = 'interBK_nC2_K4x2';
PB.nK = 4; PB.m = 2; PB.M_K = 16;
PB.lambda      = 0;
PB.sigma       = 0.14999999999999999;
PB.b           = 0.059999999999999998;
PB.q           = 0;
PB.j_min       = 0;
PB.obstacle_ratio = 0;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 3;
PB.relabel     = 1;
PB.lb_orig     = zeros(6,1); PB.ub_orig = ones(6,1);
PB.lb_C_work   = PB.lb_orig; PB.ub_C_work = PB.ub_orig;
PB.ProblemData = ProblemData;
PB.cardinality = [Inf Inf [2 2 2 2]];
PB.cardinality = PB.cardinality(:)';
PB.n_finite    = 4;
PB.finite_part_size = 16;
PB.domain_size = Inf;
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by construction
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.952392578125 0.96768785246151501 2 2 2 2];   % protocol starting point: Halton point 3023 (index representation)
PB.x0_alt_idx  = [0.438232421875 0.89772900472488948 1 2 2 1;0.696044921875 0.95945739978661781 2 2 1 1;0.192138671875 0.9772900472488949 2 2 1 1;0.824951171875 0.90733119951226948 1 2 2 2;0.315185546875 0.93613778387440927 2 2 2 2;0.572998046875 0.99786617893613772 1 2 1 2;0.069091796875 0.91693339429964948 2 2 2 2;0.889404296875 0.94573997866178938 1 2 2 1;0.376708984375 0.97454656302392928 1 2 1 1];   % 9 alternative starting points (Halton indices 3023 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K4x2_s3_b0p06_l0: n = 6 (nC = 2, nK = 4 x 2), lambda = 0, b = 0.06, j_min = 0, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 6 || ~all(isfinite(z)), error('interBK_nC2_K4x2_s3_b0p06_l0: z must be a finite numeric vector with 6 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K4x2_s3_b0p06_l0: continuous coordinates must lie in [0,1].'); end
kk = z(3:6);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 2), error('interBK_nC2_K4x2_s3_b0p06_l0: label indices must be integers in 1..2.'); end
h = 0;
for i = 1:4, h = h + (D.l2s{i}(kk(i)) == D.ks(i)); end
if h == 0, shift = D.lambda * D.s; else, shift = zeros(2,1); end
dd = u - D.xs - shift;
f = dd' * dd + D.b * (1 - h / 4);
end

function xr = idx2cell(z, PB)
xr = cell(1, 6);
for i = 1:2, xr{i} = z(i); end
for i = 3:6
    k = z(i); if k ~= floor(k) || k < 1 || k > 2, error('interBK_nC2_K4x2_s3_b0p06_l0: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 6, error('interBK_nC2_K4x2_s3_b0p06_l0: xr must be a cell array with 6 entries.'); end
z = zeros(1, 6);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K4x2_s3_b0p06_l0: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:6
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K4x2_s3_b0p06_l0: unknown label %s.', lab); end
    z(i) = j;
end
end
