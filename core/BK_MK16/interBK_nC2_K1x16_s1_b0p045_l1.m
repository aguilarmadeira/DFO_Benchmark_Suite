function PB = interBK_nC2_K1x16_s1_b0p045_l1()
%INTERBK_NC2_K1X16_S1_B0P045_L1  DFO Benchmark Suite v2 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 1 categorical variables with m = 16 categories each (M_K = 16 combinations)
% lambda = 1, sigma = 0.15, b = 0.045, q = lambda^2 sigma^2 / b = 0.5, j_min = 1 (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.66935997009277348 0.31318170784715516];
s_ref = [-0.077213565594290476;0.12860040936255296];
l2s = {[15 1 3 4 13 9 12 14 6 0 11 8 10 7 5 2]};
k_ref = 14;
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2';'k3_3';'k3_4';'k3_5';'k3_6';'k3_7';'k3_8';'k3_9';'k3_10';'k3_11';'k3_12';'k3_13';'k3_14';'k3_15';'k3_16'}};
D = struct('name', 'interBK_nC2_K1x16_s1_b0p045_l1', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 1, 'b', 0.044999999999999998);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K1x16_s1_b0p045_l1';
PB.base_function = 'familyBK';
PB.family      = 'BK';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core, Family BK (controlled categorical dimension): 1 categorical variables x 16 categories.';
PB.nvar        = 3;
PB.nobj        = 1;
PB.var_types   = {'C','C','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 1);
PB.proportion  = [2 0 1];
PB.pattern     = 'BK';
PB.configuration = 'interBK_nC2_K1x16';
PB.nK = 1; PB.m = 16; PB.M_K = 16;
PB.lambda      = 1;
PB.sigma       = 0.14999999999999999;
PB.b           = 0.044999999999999998;
PB.q           = 0.5;
PB.j_min       = 1;
PB.obstacle_ratio = 0.5;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 1;
PB.relabel     = 1;
PB.lb_orig     = zeros(3,1); PB.ub_orig = ones(3,1);
PB.lb_C_work   = PB.lb_orig; PB.ub_C_work = PB.ub_orig;
PB.ProblemData = ProblemData;
PB.cardinality = [Inf Inf 16];
PB.cardinality = PB.cardinality(:)';
PB.n_finite    = 1;
PB.finite_part_size = 16;
PB.domain_size = Inf;
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by construction
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.85888671875 0.50983081847279388 5];   % protocol starting point: Halton point 2011 (index representation)
PB.x0_alt_idx  = [0.343994140625 0.53863740283493367 16;0.601806640625 0.46867855509830819 8;0.097900390625 0.48163389727175732 1;0.918212890625 0.54336229233348576 11;0.408447265625 0.47340344459686018 4;0.666259765625 0.50221002895900013 14;0.162353515625 0.52004267642127722 7;0.795166015625 0.45008382868465169 3;0.282470703125 0.51181222374638013 10];   % 9 alternative starting points (Halton indices 2011 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K1x16_s1_b0p045_l1: n = 3 (nC = 2, nK = 1 x 16), lambda = 1, b = 0.045, j_min = 1, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('interBK_nC2_K1x16_s1_b0p045_l1: z must be a finite numeric vector with 3 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K1x16_s1_b0p045_l1: continuous coordinates must lie in [0,1].'); end
kk = z(3:3);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 16), error('interBK_nC2_K1x16_s1_b0p045_l1: label indices must be integers in 1..16.'); end
h = 0;
for i = 1:1, h = h + (D.l2s{i}(kk(i)) == D.ks(i)); end
if h == 0, shift = D.lambda * D.s; else, shift = zeros(2,1); end
dd = u - D.xs - shift;
f = dd' * dd + D.b * (1 - h / 1);
end

function xr = idx2cell(z, PB)
xr = cell(1, 3);
for i = 1:2, xr{i} = z(i); end
for i = 3:3
    k = z(i); if k ~= floor(k) || k < 1 || k > 16, error('interBK_nC2_K1x16_s1_b0p045_l1: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 3, error('interBK_nC2_K1x16_s1_b0p045_l1: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K1x16_s1_b0p045_l1: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:3
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K1x16_s1_b0p045_l1: unknown label %s.', lab); end
    z(i) = j;
end
end
