function PB = shubert_d_2D()
%SHUBERT_D_2D  DFO Benchmark Suite v2.1.0 core: shubert, dimension 2, pattern D, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -186.7309088310238 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:D
%
lb_orig = [-10 -6.3071991058273635];
ub_orig = [7.1274169921874995 10];
lb_C_work = [-10 -6.3071991058273635];
ub_C_work = [7.1274169921874995 10];
T_break = [-7.0312477213541662 4.8904109468407597];   % breakpoints of T (NaN for continuous)
T_image = [-7.0835064071999998 4.8580568773000001];   % their images
m = 151;
ProblemData = cell(1, 2);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = {linspace(lb_orig(2), ub_orig(2), m)'};
PB.name        = 'shubert_d_2D';
PB.base_function = 'shubert';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: D adaptation of shubert (F = f o T).';
PB.nvar        = 2;
PB.nobj        = 1;
PB.var_types   = {'D','D'};
PB.type_counts = struct('C', 0, 'D', 2, 'K', 0);
PB.proportion  = [0 2 0];
PB.pattern     = 'D';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [151 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 1;
PB.domain_size = 22801;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shubert_d_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN];
D.fn = @(xh) shubert_orig(xh);
PB.f_star      = -186.73090883102384;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 12;
PB.x0_idx      = [34 149];   % protocol starting point: Halton point 1052 (index representation)
PB.x0_alt_idx  = [129 139;53 143;92 146;15 135;138 145;63 149;102 136;25 141;120 150];   % 9 alternative starting points (Halton indices 1052 + 63 k)
PB.core_pidx   = 52;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shubert_d_2D: n = 2, pattern D, m = 151, f* = -186.7309088 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('shubert_d_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.lb(1) + (k - 1) * (D.ub(1) - D.lb(1)) / (D.m - 1);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
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
xr = cell(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 2, error('shubert_d_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shubert_d_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shubert_d_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shubert_d_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shubert_d_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = shubert_orig(x)
%
% Purpose:
%
%    Function shubert is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -186.7309
%    Global minima (18)
%    Local minima (760)
%    Search domain: -10 <= x <= 10 (Huyer and Neumaier (1999))
%                                  (Huyer and Neumaier (2008))
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
f = sum([1:5].*cos([2:6].*x(1)+[1:5]),2)*sum([1:5].*cos([2:6].*x(2)+[1:5]),2);
end
%
% End of shubert.

