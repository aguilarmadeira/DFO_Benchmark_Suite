function PB = interBK_nC2_K2x4_s2_b0p03_l1()
%INTERBK_NC2_K2X4_S2_B0P03_L1  DFO Benchmark Suite v2 core, Family BK: f = ||x_C - x_C* - lambda s(kappa)||^2 + b (1 - h(kappa)/nK).
%
% Copyright (C) 2026 J. F. A. Madeira
%
% nC = 2 continuous variables in [0,1]; nK = 2 categorical variables with m = 4 categories each (M_K = 16 combinations)
% lambda = 1, sigma = 0.15, b = 0.03, q = lambda^2 sigma^2 / b = 0.75, j_min = 2 (minimal number of simultaneous categorical
% changes from an antipodal conditional minimiser that improves f without continuous re-optimisation; Inf: none; 0: lambda = 0)
% h(kappa) = number of categorical coordinates equal to the optimal combination; s(kappa) = s if h = 0, else 0.
% Global optimum f* = 0 (E: by construction); the optimal solution is listed in reference_optima.csv of this
% folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
%
x_C_ref = [0.73419418334960929 0.4307744415654795];
s_ref = [-0.14980782553367489;-0.0075904814650998775];
l2s = {[3 1 2 0], [1 3 0 2]};
k_ref = [2 3];
ProblemData = {[0; 1], [0; 1], {'k3_1';'k3_2';'k3_3';'k3_4'}, {'k4_1';'k4_2';'k4_3';'k4_4'}};
D = struct('name', 'interBK_nC2_K2x4_s2_b0p03_l1', 'xs', x_C_ref(:), 's', s_ref(:), 'ks', k_ref, 'lambda', 1, 'b', 0.029999999999999999);
D.l2s = l2s;
PB.name        = 'interBK_nC2_K2x4_s2_b0p03_l1';
PB.base_function = 'familyBK';
PB.family      = 'BK';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core, Family BK (controlled categorical dimension): 2 categorical variables x 4 categories.';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','C','K','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'BK';
PB.configuration = 'interBK_nC2_K2x4';
PB.nK = 2; PB.m = 4; PB.M_K = 16;
PB.lambda      = 1;
PB.sigma       = 0.14999999999999999;
PB.b           = 0.029999999999999999;
PB.q           = 0.75;
PB.j_min       = 2;
PB.obstacle_ratio = 0.75;   % f(all nK changes at fixed x_C) / f(antipodal conditional minimum) = q
PB.seed        = 2;
PB.relabel     = 1;
PB.lb_orig     = zeros(4,1); PB.ub_orig = ones(4,1);
PB.lb_C_work   = PB.lb_orig; PB.ub_C_work = PB.ub_orig;
PB.ProblemData = ProblemData;
PB.cardinality = [Inf Inf [4 4]];
PB.cardinality = PB.cardinality(:)';
PB.n_finite    = 2;
PB.finite_part_size = 16;
PB.domain_size = Inf;
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by construction
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.40576171875 0.32464563328760859 3 4];   % protocol starting point: Halton point 2022 (index representation)
PB.x0_alt_idx  = [0.640869140625 0.25468678555098306 1 4;0.148681640625 0.28349336991312302 3 4;0.769775390625 0.29644871208657214 2 4;0.277587890625 0.22648986434994664 4 4;0.517822265625 0.288218259411675 3 4;0.025634765625 0.31702484377381496 1 4;0.974853515625 0.23609205913732662 3 4;0.482666015625 0.26489864349946651 2 4;0.719970703125 0.3266270385611949 4 4];   % 9 alternative starting points (Halton indices 2022 + 63 k)
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
fprintf('interBK_nC2_K2x4_s2_b0p03_l1: n = 4 (nC = 2, nK = 2 x 4), lambda = 1, b = 0.03, j_min = 2, f* = 0 (E)\n');
end

function f = objective(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('interBK_nC2_K2x4_s2_b0p03_l1: z must be a finite numeric vector with 4 entries.'); end
u = z(1:2); u = u(:);
if any(u < 0) || any(u > 1), error('interBK_nC2_K2x4_s2_b0p03_l1: continuous coordinates must lie in [0,1].'); end
kk = z(3:4);
if any(kk ~= floor(kk)) || any(kk < 1) || any(kk > 4), error('interBK_nC2_K2x4_s2_b0p03_l1: label indices must be integers in 1..4.'); end
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
    k = z(i); if k ~= floor(k) || k < 1 || k > 4, error('interBK_nC2_K2x4_s2_b0p03_l1: label index out of range.'); end
    xr{i} = PB.ProblemData{i}{k};
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 4, error('interBK_nC2_K2x4_s2_b0p03_l1: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:2
    v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v < 0 || v > 1, error('interBK_nC2_K2x4_s2_b0p03_l1: x%d must be a finite scalar in [0,1].', i); end
    z(i) = double(v);
end
for i = 3:4
    lab = xr{i}; if isstring(lab), lab = char(lab); end
    j = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(j), error('interBK_nC2_K2x4_s2_b0p03_l1: unknown label %s.', lab); end
    z(i) = j;
end
end
