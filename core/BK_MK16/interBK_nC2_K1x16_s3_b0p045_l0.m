function PB = interBK_nC2_K1x16_s3_b0p045_l0()
%INTERBK_NC2_K1X16_S3_B0P045_L0  DFO Benchmark Suite v2.1.0 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 1 categorical variables with m = 16 categories each (M_K = 16 combinations)
% lambda = 0, sigma = 0.15, b = 0.045, q = lambda^2 sigma^2 / b = 0, j_min = 0 (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.61149864196777348 0.69269702563407765];
s_ref = [0.098328318961646663;-0.11327639511467813];
l2s = {[4 15 10 5 9 11 3 14 8 7 0 12 13 1 2 6]};
k_ref = 5;
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2';'k3_3';'k3_4';'k3_5';'k3_6';'k3_7';'k3_8';'k3_9';'k3_10';'k3_11';'k3_12';'k3_13';'k3_14';'k3_15';'k3_16'}};
D = struct('name', 'interBK_nC2_K1x16_s3_b0p045_l0', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 0, 'b', 0.044999999999999998);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K1x16_s3_b0p045_l0';
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
PB.lambda      = 0;
PB.sigma       = 0.14999999999999999;
PB.b           = 0.044999999999999998;
PB.q           = 0;
PB.j_min       = 0;
PB.obstacle_ratio = 0;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 3;
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
PB.x0_idx      = [0.73388671875 0.28760859625057156 11];   % protocol starting point: Halton point 2013 (index representation)
PB.x0_alt_idx  = [0.218994140625 0.31641518061271146 4;0.851806640625 0.24645633287608593 14;0.347900390625 0.25941167504953511 7;0.605712890625 0.3211400701112635 3;0.095947265625 0.25118122237463797 10;0.916259765625 0.27998780673677792 6;0.412353515625 0.29782045419905501 13;0.670166015625 0.22786160646242951 9;0.157470703125 0.28959000152415787 2];   % 9 alternative starting points (Halton indices 2013 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K1x16_s3_b0p045_l0: n = 3 (nC = 2, nK = 1 x 16), lambda = 0, b = 0.045, j_min = 0, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('interBK_nC2_K1x16_s3_b0p045_l0: z must be a finite numeric vector with 3 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K1x16_s3_b0p045_l0: continuous coordinates must lie in [0,1].'); end
kk = z(3:3);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 16), error('interBK_nC2_K1x16_s3_b0p045_l0: label indices must be integers in 1..16.'); end
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
    k = z(i); if k ~= floor(k) || k < 1 || k > 16, error('interBK_nC2_K1x16_s3_b0p045_l0: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 3, error('interBK_nC2_K1x16_s3_b0p045_l0: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K1x16_s3_b0p045_l0: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:3
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K1x16_s3_b0p045_l0: unknown label %s.', lab); end
    z(i) = j;
end
end
