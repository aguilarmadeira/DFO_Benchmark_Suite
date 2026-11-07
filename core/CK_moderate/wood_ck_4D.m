function PB = wood_ck_4D()
%WOOD_CK_4D  DFO Benchmark Suite v2.1.0 core: wood, dimension 4, pattern CK, strategy moderate.
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
lb_orig = [-10 -10 -10 -10];
ub_orig = [6.1211669921874998 6.5664583650866231 6.4227327999999995 7.3514369012911285];
lb_C_work = [-10 -10 -1000 -10];
ub_C_work = [6.1211669921875007 6.5664583650866231 642.27327999999989 7.3514369012911285];
T_break = [NaN 1.0443055767244172 NaN 0.98924337081771618];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 NaN 1];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'wood_ck_4D';
PB.base_function = 'wood';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CK adaptation of wood (F = f o T).';
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
D = struct('name', 'wood_ck_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [1 NaN 100 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [1 NaN 1 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([30 135 36 60 80 141 142 61 111 151 55 102 78 62 148 43 45 70 53 116 33 103 34 118 100 42 82 146 128 17 115 139 39 73 86 88 108 68 2 83 46 145 26 98 105 147 64 28 38 87 110 91 72 127 48 94 93 112 51 15 130 3 10 89 52 150 25 49 92 96 41 16 132 84 47 6 138 131 27 1 44 74 124 4 54 125 50 134 24 90 31 58 123 8 37 104 109 114 20 59 79 63 5 149 121 101 56 75 137 57 22 40 76 35 106 117 85 12 32 119 81 65 18 144 136 29 99 13 19 122 67 107 23 71 95 11 120 140 77 126 21 129 7 143 133 66 14 69 97 9 113])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([98 25 82 11 150 140 59 123 22 62 41 118 37 100 90 88 66 148 86 113 114 5 81 110 151 138 39 108 14 129 131 48 104 84 63 103 137 78 4 74 46 24 64 111 127 94 93 101 96 2 126 124 116 144 75 76 21 42 112 18 105 80 83 27 55 146 130 52 61 35 142 20 85 8 134 135 34 71 32 6 49 79 139 149 60 147 56 77 30 106 125 65 26 132 109 97 40 128 29 68 122 73 15 3 92 115 133 87 45 70 36 54 91 50 7 117 141 145 51 19 102 89 43 33 38 57 143 44 53 28 107 16 136 119 10 67 47 58 13 69 121 95 1 23 31 12 17 99 9 72 120])';
D.fn = @(xh) wood_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [4.3657860159873962 57 151.43064207359998 142];   % protocol starting point: Halton point 1063 (index representation)
PB.x0_alt_idx  = [-3.5688508629798887 67 -660.5092675584001 149;0.39846757650375331 56 482.5129353215998 134;-7.5361693024635308 59 -342.03963310080007 141;2.3821267962455757 63 -775.59977902080004 147;-5.5525100827217102 53 -10.957339852799919 130;-1.5851916432380673 62 -444.51748577280011 136;-9.5198285222053531 64 320.12495339520001 143;5.6409955143928538 54 -113.43519252480007 150;-2.293641364574432 58 -925.37510215680004 135];   % 9 alternative starting points (Halton indices 1063 + 63 k)
PB.core_pidx   = 63;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('wood_ck_4D: n = 4, pattern CK, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('wood_ck_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('wood_ck_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('wood_ck_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
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
if ~iscell(xr) || numel(xr) ~= 4, error('wood_ck_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('wood_ck_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('wood_ck_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('wood_ck_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('wood_ck_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = wood_orig(x)
%
% Purpose:
%
%    Function wood is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [1 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = 100*(x(1,1)^2 - x(2,1))^2 + (1-x(1,1))^2 + 90*(x(3,1)^2-x(4,1))^2 +...
   (1-x(3,1))^2 + 10.1*((1-x(2,1))^2 + (1-x(4,1))^2) + 19.8*(1-x(2,1))*(1-x(4,1));
end
%
% End of wood.

