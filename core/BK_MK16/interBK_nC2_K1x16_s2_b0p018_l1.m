function PB = interBK_nC2_K1x16_s2_b0p018_l1()
%INTERBK_NC2_K1X16_S2_B0P018_L1  DFO Benchmark Suite v2.1.0 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 1 categorical variables with m = 16 categories each (M_K = 16 combinations)
% lambda = 1, sigma = 0.15, b = 0.018, q = lambda^2 sigma^2 / b = 1.25, j_min = Inf (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.27590293884277345 0.45000000000000001];
s_ref = [0.07188053468826594;-0.13165556856027394];
l2s = {[8 0 6 2 1 15 5 9 7 13 12 4 14 3 11 10]};
k_ref = 9;
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2';'k3_3';'k3_4';'k3_5';'k3_6';'k3_7';'k3_8';'k3_9';'k3_10';'k3_11';'k3_12';'k3_13';'k3_14';'k3_15';'k3_16'}};
D = struct('name', 'interBK_nC2_K1x16_s2_b0p018_l1', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 1, 'b', 0.017999999999999999);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K1x16_s2_b0p018_l1';
PB.base_function = 'familyBK';
PB.family      = 'BK';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core, Family BK (controlled categorical dimension): 1 categorical variables x 16 categories.';
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
PB.b           = 0.017999999999999999;
PB.q           = 1.25;
PB.j_min       = Inf;
PB.obstacle_ratio = 1.25;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 2;
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
PB.x0_idx      = [0.23388671875 0.84316415180612703 8];   % protocol starting point: Halton point 2012 (index representation)
PB.x0_alt_idx  = [0.843994140625 0.87197073616826692 1;0.351806640625 0.80201188843164151 11;0.597900390625 0.81496723060509058 4;0.105712890625 0.87669562566681902 14;0.908447265625 0.80673677793019349 7;0.416259765625 0.83554336229233339 3;0.662353515625 0.85337600975461048 10;0.170166015625 0.78341716201798495 6;0.782470703125 0.84514555707971328 13];   % 9 alternative starting points (Halton indices 2012 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K1x16_s2_b0p018_l1: n = 3 (nC = 2, nK = 1 x 16), lambda = 1, b = 0.018, j_min = Inf, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('interBK_nC2_K1x16_s2_b0p018_l1: z must be a finite numeric vector with 3 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K1x16_s2_b0p018_l1: continuous coordinates must lie in [0,1].'); end
kk = z(3:3);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 16), error('interBK_nC2_K1x16_s2_b0p018_l1: label indices must be integers in 1..16.'); end
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
    k = z(i); if k ~= floor(k) || k < 1 || k > 16, error('interBK_nC2_K1x16_s2_b0p018_l1: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 3, error('interBK_nC2_K1x16_s2_b0p018_l1: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K1x16_s2_b0p018_l1: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:3
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K1x16_s2_b0p018_l1: unknown label %s.', lab); end
    z(i) = j;
end
end
