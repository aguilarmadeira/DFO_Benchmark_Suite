function PB = cosine_mixture_cdk_4D()
%COSINE_MIXTURE_CDK_4D  DFO Benchmark Suite v2 core: cosine_mixture, dimension 4, pattern CDK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -4.4 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C
%
lb_orig = [-1 -1 -1 -1];
ub_orig = [1 1 1 1];
lb_C_work = [-1 -1 -1 -100];
ub_C_work = [1 1 1 100];
T_break = [NaN 1 1 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 1 NaN];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
PB.name        = 'cosine_mixture_cdk_4D';
PB.base_function = 'cosine_mixture';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CDK adaptation of cosine_mixture (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C'};
PB.type_counts = struct('C', 2, 'D', 1, 'K', 1);
PB.proportion  = [2 1 1];
PB.pattern     = 'CDK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100];
PB.cardinality = [Inf 151 151 Inf];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 151;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'cosine_mixture_cdk_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([38 121 109 75 118 138 137 65 148 147 56 44 35 108 62 26 77 93 17 60 55 33 3 94 36 90 53 98 66 23 103 101 115 119 6 8 124 145 122 149 82 80 37 68 2 150 16 54 141 72 96 81 15 89 128 144 88 70 7 19 39 105 50 58 112 21 29 84 135 120 79 92 86 146 12 125 43 110 116 113 57 74 130 42 69 52 73 85 18 123 102 97 27 71 134 142 91 48 139 106 31 76 127 46 107 117 131 59 95 63 104 9 34 20 64 129 67 49 4 1 51 22 11 13 30 140 100 136 151 83 99 5 132 87 126 32 61 14 111 47 24 40 45 41 114 78 143 28 25 10 133])';
D.fn = @(xh) cosine_mixture_orig(xh);
PB.f_star      = -4.4000000000000004;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 16;
PB.x0_idx      = [0.123046875 59 128 -51.436901291128699];   % protocol starting point: Halton point 1009 (index representation)
PB.x0_alt_idx  = [-0.9052734375 63 88 -46.68887963348606;0.9228515625 53 14 -65.930862140774678;-0.0849609375 62 118 -57.184506455643486;0.4306640625 64 43 -48.438150770512287;-0.5888671875 54 147 -67.68013327780092;0.6767578125 58 73 -62.932111620158274;-0.3310546875 68 5 -54.185755935027082;0.1845703125 55 104 -45.43940024989589;-0.8408203125 59 35 -64.681382757184508];   % 9 alternative starting points (Halton indices 1009 + 63 k)
PB.core_pidx   = 9;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('cosine_mixture_cdk_4D: n = 4, pattern CDK, m = 151, f* = -4.4 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('cosine_mixture_cdk_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('cosine_mixture_cdk_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('cosine_mixture_cdk_4D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
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
xr = cell(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 4, error('cosine_mixture_cdk_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('cosine_mixture_cdk_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('cosine_mixture_cdk_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('cosine_mixture_cdk_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('cosine_mixture_cdk_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = cosine_mixture_orig(x)
%
% Purpose:
%
%    Function cosine_mixture is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value: -2 if n = 2
%                          -4 if n = 4
%    Global minima
%    Local minima
%    Search domain: -1 <= x <= 1 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
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
f = 0.1*sum(cos(5*pi.*x),1) - sum(x.^2,1); 
end
%
% End of cosine_mixture.

