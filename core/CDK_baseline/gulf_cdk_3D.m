function PB = gulf_cdk_3D()
%GULF_CDK_3D  DFO Benchmark Suite v2.1.0 core: gulf, dimension 3, pattern CDK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 9.709124424797534e-31 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D
%
lb_orig = [0.10000000000000001 4.1859895341157349 0];
ub_orig = [86.607522094726562 25.600000000000001 4.0132832000000001];
lb_C_work = [0.10000000000000001 4.1859895341157349 0];
ub_C_work = [86.607522094726562 25.600000000000001 4.0132832000000001];
T_break = [NaN 25.028959720909757 1.4982923946666666];   % breakpoints of T (NaN for continuous)
T_image = [NaN 25 1.5];   % their images
m = 151;
ProblemData = cell(1, 3);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
PB.name        = 'gulf_cdk_3D';
PB.base_function = 'gulf';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CDK adaptation of gulf (F = f o T).';
PB.nvar        = 3;
PB.nobj        = 1;
PB.var_types   = {'C','K','D'};
PB.type_counts = struct('C', 1, 'D', 1, 'K', 1);
PB.proportion  = [1 1 1];
PB.pattern     = 'CDK';
PB.strategy    = 'baseline';
PB.scale_factors_C = 1;
PB.cardinality = [Inf 151 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 151;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'gulf_cdk_3D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 3);
D.cw = [50.000000000000007 NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [50 NaN NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([114 37 71 63 52 120 111 19 119 29 30 113 51 32 90 129 43 26 118 27 20 1 49 123 58 80 44 22 146 98 96 47 9 13 83 86 48 61 36 66 82 75 115 85 24 121 41 45 109 102 23 128 116 97 151 35 144 57 64 21 125 112 60 130 73 103 140 81 4 33 74 138 25 117 142 92 133 105 3 147 110 5 16 134 18 141 42 101 84 8 89 15 126 12 31 131 59 28 50 76 56 108 40 106 88 127 34 10 53 67 122 145 148 11 136 135 139 104 78 14 107 7 137 143 46 99 100 77 72 149 91 70 95 69 79 94 54 132 39 93 17 68 87 124 55 62 150 2 6 38 65])';
D.fn = @(xh) gulf_orig(xh);
PB.f_star      = 9.7091244247975338e-31;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.1422400010228157 98 146];   % protocol starting point: Halton point 1024 (index representation)
PB.x0_alt_idx  = [85.29808206301928 87 77;42.720161032021046 91 3;64.009121547520166 101 108;21.431200516521933 89 32;74.653601805269716 93 137;32.07568077427149 97 63;53.36464128977061 84 23;10.786720258772373 94 93;79.975841934144498 98 53];   % 9 alternative starting points (Halton indices 1024 + 63 k)
PB.core_pidx   = 24;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('gulf_cdk_3D: n = 3, pattern CDK, m = 151, f* = 9.709124425e-31 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('gulf_cdk_3D: z must be a finite numeric vector with 3 entries.'); end
xh = zeros(3, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('gulf_cdk_3D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
end

function x = T_map(g, a, b, zs, xs)
% piecewise-affine T with T(zs) = xs; fixes the bounds when zs is interior
if abs(g - zs) <= 1e-12 * max(1, abs(b - a)), x = xs; return; end
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
xr = cell(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 3, error('gulf_cdk_3D: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('gulf_cdk_3D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('gulf_cdk_3D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('gulf_cdk_3D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('gulf_cdk_3D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = gulf_orig(x)
%
% Purpose:
%
%    Function gulf is the function described in
%    Montaz et al. (2005).
%
%    dim = 3
%    Minimum global value = 0
%    Global minima = [50 25 1.5]'
%    Local minima 
%    Search domain: 0.1 <= x(1) <= 100 
%                   0 <= x(2) <= 25.6
%                   0 <= x(3) <= 5    (Montaz et al. (2005))
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
u = 25+(-50*log(0.01.*[1:99])).^(2/3);
f = sum((exp(-(u-x(2)).^x(3)./x(1))-0.01.*[1:99]).^2,2);
end
%
% End of gulf.

