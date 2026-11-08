function PB = schaffer2_ck_2D()
%SCHAFFER2_CK_2D  DFO Benchmark Suite v2.1.0 core: schaffer2, dimension 2, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 65 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K
%
lb_orig = [-100 -64.704364172128237];
ub_orig = [65.149169921875 100];
lb_C_work = [-100 -64.704364172128237];
ub_C_work = [65.149169921875 100];
T_break = [NaN -0.36672191739064885];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0];   % their images
m = 65;
ProblemData = cell(1, 2);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'schaffer2_ck_2D';
PB.base_function = 'schaffer2';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CK adaptation of schaffer2 (F = f o T).';
PB.nvar        = 2;
PB.nobj        = 1;
PB.var_types   = {'C','K'};
PB.type_counts = struct('C', 1, 'D', 0, 'K', 1);
PB.proportion  = [1 0 1];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = 1;
PB.cardinality = [Inf 65];
PB.n_finite    = 1;
PB.finite_part_size = 65;
PB.M_K = 65;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'schaffer2_ck_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.cw = [0 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [0 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([60 24 32 25 63 36 19 65 9 54 45 22 59 42 44 53 61 2 57 29 14 15 58 21 12 64 30 17 34 56 1 35 18 23 41 3 48 28 49 51 33 20 47 37 27 11 5 43 50 62 38 16 52 40 46 39 10 4 13 8 26 6 31 7 55])';
D.fn = @(xh) schaffer2_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [8.4597820043563843 29];   % protocol starting point: Halton point 1045 (index representation)
PB.x0_alt_idx  = [-72.824575066566467 24;29.748542189598083 26;-51.535814881324768 27;-11.861307263374329 23;-93.14566433429718 27;55.875656962394714 29;-25.408700108528137 23;13.781972050666809 25;-67.502385020256042 29];   % 9 alternative starting points (Halton indices 1045 + 63 k)
PB.core_pidx   = 45;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('schaffer2_ck_2D: n = 2, pattern CK, m = 65, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('schaffer2_ck_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('schaffer2_ck_2D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
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
if ~iscell(xr) || numel(xr) ~= 2, error('schaffer2_ck_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('schaffer2_ck_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('schaffer2_ck_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('schaffer2_ck_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('schaffer2_ck_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = schaffer2_orig(x)
%
% Purpose:
%
%    Function schaffer2 is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -100 <= x <= 100 (Montaz et al. (2005))
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
f = (x(1)^2+x(2)^2)^0.25*(sin(50*(x(1)^2+x(2)^2)^0.1)^2+1);
end
%
% End of schaffer2.

