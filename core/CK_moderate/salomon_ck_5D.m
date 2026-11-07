function PB = salomon_ck_5D()
%SALOMON_CK_5D  DFO Benchmark Suite v2 core: salomon, dimension 5, pattern CK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K x5:C
%
lb_orig = [-100 -100 -100 -100 -100];
ub_orig = [100 100 100 100 100];
lb_C_work = [-100 -100 -10000 -100 -1];
ub_C_work = [100 100 10000 100 1];
T_break = [NaN 0 NaN 0 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 NaN 0 NaN];   % their images
m = 151;
ProblemData = cell(1, 5);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
PB.name        = 'salomon_ck_5D';
PB.base_function = 'salomon';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of salomon (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C'};
PB.type_counts = struct('C', 3, 'D', 0, 'K', 2);
PB.proportion  = [3 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01];
PB.cardinality = [Inf 151 Inf 151 Inf];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'salomon_ck_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([116 20 59 151 101 142 94 71 88 72 79 107 114 52 45 106 73 99 124 38 150 10 105 61 113 32 89 126 111 37 98 26 64 2 42 7 16 80 75 34 122 82 67 115 23 66 44 55 109 112 33 143 21 135 119 90 14 93 125 147 30 85 47 149 132 36 49 129 103 43 51 104 50 68 146 96 1 31 131 128 63 87 18 84 41 60 127 22 56 5 145 35 95 54 148 123 39 24 133 117 4 137 8 83 12 121 130 15 69 29 91 65 118 92 144 28 6 139 48 53 76 141 57 3 40 9 100 27 138 19 81 134 102 17 25 86 74 77 140 120 62 97 70 108 13 110 136 78 58 46 11])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([141 77 104 69 15 99 40 145 103 48 89 128 96 36 106 88 79 26 60 140 92 38 13 30 74 41 136 59 101 70 25 114 112 109 7 9 44 144 107 118 137 130 49 17 53 126 46 127 115 19 78 116 117 124 51 8 100 34 76 3 135 43 50 94 91 39 52 4 82 122 55 22 54 151 98 62 102 35 119 32 149 90 129 5 121 66 133 29 142 11 20 110 95 61 86 134 47 33 83 16 14 84 18 93 23 148 31 27 123 42 105 80 68 108 138 45 24 71 125 73 97 85 81 28 146 37 2 132 111 147 120 58 12 65 75 6 1 10 56 143 21 63 150 131 57 67 139 87 113 64 72])';
D.fn = @(xh) salomon_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-43.65234375 94 662.40000000000146 133 0.56574004507888831];   % protocol starting point: Halton point 1042 (index representation)
PB.x0_alt_idx  = [7.91015625 98 -8457.6000000000004 139 -0.060856498873027798;-92.87109375 88 4694.4000000000015 146 -0.50713749060856506;89.94140625 90 -4425.5999999999995 132 0.84973703981968463;-12.01171875 100 8572.7999999999993 138 0.40345604808414737;39.55078125 89 -547.20000000000073 145 -0.23953281879653032;-61.23046875 93 -6595.1999999999998 148 -0.68581381053206747;65.33203125 96 3484.8000000000029 134 0.67106071989618221;-37.20703125 85 -2563.2000000000007 140 0.22477972816064473;14.35546875 94 7516.8000000000029 147 -0.40181681579127115];   % 9 alternative starting points (Halton indices 1042 + 63 k)
PB.core_pidx   = 42;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('salomon_ck_5D: n = 5, pattern CK, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('salomon_ck_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('salomon_ck_5D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('salomon_ck_5D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('salomon_ck_5D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
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
xr = cell(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 5, error('salomon_ck_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('salomon_ck_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('salomon_ck_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('salomon_ck_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('salomon_ck_5D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = salomon_orig(x)
%
% Purpose:
%
%    Function salomon is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
%    Cases considered: n = 5,10 Montaz et al. (2005)
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
f = 1-cos(2*pi*norm(x,2))+0.1*norm(x,2);
end
%
% End of salomon.

