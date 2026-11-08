function PB = gulf_d_3D()
%GULF_D_3D  DFO Benchmark Suite v2.1.0 core: gulf, dimension 3, pattern D, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 9.709124424797534e-31 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:D x2:D x3:D
%
lb_orig = [0.10000000000000001 4.1859895341157349 0];
ub_orig = [86.607522094726562 25.600000000000001 4.0132832000000001];
lb_C_work = [0.10000000000000001 4.1859895341157349 0];
ub_C_work = [86.607522094726562 25.600000000000001 4.0132832000000001];
T_break = [49.841825204467781 25.064649738352895 1.5049812];   % breakpoints of T (NaN for continuous)
T_image = [50 25 1.5];   % their images
m = 41;
ProblemData = cell(1, 3);
ProblemData{1} = {linspace(lb_orig(1), ub_orig(1), m)'};
ProblemData{2} = {linspace(lb_orig(2), ub_orig(2), m)'};
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
PB.name        = 'gulf_d_3D';
PB.base_function = 'gulf';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: D adaptation of gulf (F = f o T).';
PB.nvar        = 3;
PB.nobj        = 1;
PB.var_types   = {'D','D','D'};
PB.type_counts = struct('C', 0, 'D', 3, 'K', 0);
PB.proportion  = [0 3 0];
PB.pattern     = 'D';
PB.strategy    = 'baseline';
PB.scale_factors_C = [];
PB.cardinality = [41 41 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 1;
PB.domain_size = 68921;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'gulf_d_3D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 3);
D.cw = [NaN NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [NaN NaN NaN];
D.fn = @(xh) gulf_orig(xh);
PB.f_star      = 9.7091244247975338e-31;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [1 27 40];   % protocol starting point: Halton point 1024 (index representation)
PB.x0_alt_idx  = [40 24 21;21 25 2;31 28 29;11 24 9;35 25 37;16 27 17;26 23 7;6 26 25;38 27 15];   % 9 alternative starting points (Halton indices 1024 + 63 k)
PB.core_pidx   = 24;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('gulf_d_3D: n = 3, pattern D, m = 41, f* = 9.709124425e-31 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('gulf_d_3D: z must be a finite numeric vector with 3 entries.'); end
xh = zeros(3, 1);
k = check_idx(z(1), D.m, D.name, 1); g = D.lb(1) + (k - 1) * (D.ub(1) - D.lb(1)) / (D.m - 1);
xh(1) = T_map(g, D.lb(1), D.ub(1), D.Tb(1), D.Ti(1));
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
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
if ~iscell(xr) || numel(xr) ~= 3, error('gulf_d_3D: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('gulf_d_3D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('gulf_d_3D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('gulf_d_3D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('gulf_d_3D: x%d unknown label %s.', i, lab); end
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

