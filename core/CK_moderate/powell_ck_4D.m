function PB = powell_ck_4D()
%POWELL_CK_4D  DFO Benchmark Suite v2.1.0 core: powell, dimension 4, pattern CK, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K
%
lb_orig = [-6.0063232421875004 -6.7815475283239337 -10 -10];
ub_orig = [10 10 6.1539327999999998 6.4698042482299041];
lb_C_work = [-6.0063232421874995 -6.7815475283239337 -1000 -10];
ub_C_work = [10 10 615.39328 6.4698042482299041];
T_break = [NaN 0.042948466527799489 NaN -0.0083187560738569033];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 NaN 0];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'powell_ck_4D';
PB.base_function = 'powell';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CK adaptation of powell (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100];
PB.cardinality = [Inf 151 Inf 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'powell_ck_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [8.8817841970012523e-16 NaN 0 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [0 NaN 0 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([86 57 13 44 125 103 104 143 63 119 66 118 111 107 126 98 112 151 93 144 94 56 69 145 110 73 135 35 89 122 115 102 16 137 25 139 36 149 84 127 150 88 92 38 1 87 20 130 113 121 59 15 45 49 17 29 116 55 90 108 53 71 7 81 40 133 132 23 30 65 2 39 95 14 147 76 18 48 79 43 105 136 32 142 52 67 117 28 61 5 140 31 74 96 70 134 26 8 54 85 80 146 131 51 77 128 75 68 11 72 22 101 83 3 10 124 129 42 33 27 41 114 97 37 62 24 100 120 12 60 82 4 58 47 34 19 50 109 106 138 6 91 123 64 148 46 78 141 9 21 99])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([70 112 50 9 72 129 66 148 104 36 141 136 49 138 140 59 133 108 123 33 20 86 30 100 44 58 149 111 151 76 11 13 74 53 118 28 83 150 101 16 77 143 32 103 127 14 22 25 69 65 7 135 120 126 19 31 145 110 125 115 134 71 10 21 43 121 8 27 105 68 79 139 93 107 95 4 23 119 5 48 1 109 37 98 142 73 91 113 67 12 137 39 81 116 78 55 34 17 124 94 90 56 92 61 88 38 144 147 128 85 35 146 130 18 26 54 29 60 52 47 62 15 131 63 89 51 2 24 46 97 122 45 102 96 84 75 114 41 42 82 106 57 117 87 40 3 80 6 64 99 132])';
D.fn = @(xh) powell_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [9.0074203848838792 77 442.74004623359997 69];   % protocol starting point: Halton point 1039 (index representation)
PB.x0_alt_idx  = [1.1293081641197205 82 -293.87928944639987 75;5.0683642745018007 71 -782.37421731840004 82;-2.8097479462623594 74 31.783995801600213 67;7.0378923296928404 83 -469.11715246080007 74;-0.84021989107131922 72 345.04106065920018 80;3.098836219310761 77 -143.45386721279988 84;-4.7792760014533995 79 -880.07320289279994 69;8.0226563572883585 68 182.20941803520009 76;0.14454413652420062 78 -554.40991764480009 83];   % 9 alternative starting points (Halton indices 1039 + 63 k)
PB.core_pidx   = 39;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('powell_ck_4D: n = 4, pattern CK, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('powell_ck_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('powell_ck_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('powell_ck_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
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
if ~iscell(xr) || numel(xr) ~= 4, error('powell_ck_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('powell_ck_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('powell_ck_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('powell_ck_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('powell_ck_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = powell_orig(x)
%
% Purpose:
%
%    Function powell is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = zeros(4,1)
%    Local minima (equal to global minima)
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
f = (x(1,1) + 10*x(2,1))^2 + 5*(x(3,1) - x(4,1))^2 + (x(2,1) -...
    2*x(3,1))^4 + 10*(x(1,1) - x(4,1))^4;
end
%
% End of powell.

