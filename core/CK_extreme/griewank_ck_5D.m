function PB = griewank_ck_5D()
%GRIEWANK_CK_5D  DFO Benchmark Suite v2 core: griewank, dimension 5, pattern CK, strategy extreme.
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
lb_orig = [-600 -600 -600 -600 -600];
ub_orig = [600 600 600 600 600];
lb_C_work = [-600 -600 -60000000000 -600 -60000000000];
ub_C_work = [600 600 60000000000 600 60000000000];
T_break = [NaN 0 NaN 0 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 NaN 0 NaN];   % their images
m = 151;
ProblemData = cell(1, 5);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
PB.name        = 'griewank_ck_5D';
PB.base_function = 'griewank';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: CK adaptation of griewank (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K','C'};
PB.type_counts = struct('C', 3, 'D', 0, 'K', 2);
PB.proportion  = [3 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 100000000 100000000];
PB.cardinality = [Inf 151 Inf 151 Inf];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'griewank_ck_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([43 32 102 37 41 46 15 90 105 60 98 130 125 79 66 73 26 83 13 70 38 74 71 62 6 147 51 93 21 64 114 1 4 8 82 12 19 28 30 59 39 3 116 77 136 27 34 36 84 2 112 97 108 129 151 122 31 81 86 110 33 42 47 89 103 25 85 52 149 87 115 53 7 20 119 94 11 135 61 17 120 128 146 44 78 76 16 29 138 134 111 75 142 72 92 56 118 109 10 133 121 69 124 150 40 131 88 63 144 9 58 126 49 57 23 140 80 132 113 137 14 24 68 99 54 104 100 65 50 106 48 148 35 123 141 117 45 96 127 5 139 18 55 107 95 145 67 22 91 143 101])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([68 101 54 69 148 33 103 27 37 111 58 107 32 140 98 8 96 78 146 105 62 116 67 3 43 137 15 141 94 74 1 151 135 9 114 119 86 145 31 47 118 95 7 106 18 102 88 134 46 4 77 24 143 90 17 71 34 122 147 144 36 45 83 59 14 150 39 55 61 20 26 57 11 123 53 108 16 115 66 131 129 70 81 50 84 92 110 121 132 75 104 82 25 125 142 113 30 117 19 44 23 133 112 120 51 60 100 76 10 22 72 6 64 73 40 21 128 41 13 65 35 48 93 124 80 5 149 97 85 12 87 139 79 136 42 38 109 29 63 89 99 126 49 130 138 91 56 52 2 28 127])';
D.fn = @(xh) griewank_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [598.828125 48 31814400000.000015 44 -54320060105.184074];   % protocol starting point: Halton point 1023 (index representation)
PB.x0_alt_idx  = [-18.1640625 37 -22905599999.999992 26 37911344853.493622;291.2109375 41 56006400000 33 315552216.37866211;-313.4765625 51 1286399999.9999924 39 -37280240420.736282;445.8984375 39 -58963200000 25 54951164537.941406;-165.8203125 43 24556800000 31 16371832525.100739;143.5546875 47 -34771200000 35 -10405026979.031487;-461.1328125 34 48748800000 41 -48000819616.146439;523.2421875 44 -10579200000 27 44230585342.53125;-91.9921875 48 -42259200000 34 6634792705.4162827];   % 9 alternative starting points (Halton indices 1023 + 63 k)
PB.core_pidx   = 23;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('griewank_ck_5D: n = 5, pattern CK, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('griewank_ck_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('griewank_ck_5D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('griewank_ck_5D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
k = check_idx(z(4), D.m, D.name, 4); g = D.E{4}(k);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('griewank_ck_5D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
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
if ~iscell(xr) || numel(xr) ~= 5, error('griewank_ck_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_ck_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('griewank_ck_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('griewank_ck_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('griewank_ck_5D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = griewank_orig(x)
%
% Purpose:
%
%    Function griewank is the function described in
%    Storn and Price (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = zeros(n,1)
%    Local minima (several)
%    Search domain: -600 <= x <= 600 (Huyer and Neumaier (1999))
%                   -400 <= x <= 400 (Storn and Price (1997))
%    Cases considered: n = 10 Storn and Price (1997)
%                      n = 5 Huyer and Neumaier (1999), ICEO
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
n = size(x,1);
f = 1 + 1/4000*sum(x.^2,1)- prod(cos(x./sqrt([1:n]')),1);
end
%
% End of griewank.

