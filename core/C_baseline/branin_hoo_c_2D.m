function PB = branin_hoo_c_2D()
%BRANIN_HOO_C_2D  DFO Benchmark Suite v2 core: branin_hoo, dimension 2, pattern C, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0.3978873577297382 (A: analytic reduction, 1-D root to machine precision); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C
%
lb_orig = [-5 0];
ub_orig = [10 15];
lb_C_work = [-5 0];
ub_C_work = [10 15];
T_break = [NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 2);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
PB.name        = 'branin_hoo_c_2D';
PB.base_function = 'branin_hoo';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of branin_hoo (F = f o T).';
PB.nvar        = 2;
PB.nobj        = 1;
PB.var_types   = {'C','C'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 0);
PB.proportion  = [2 0 0];
PB.pattern     = 'C';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1];
PB.cardinality = [Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'branin_hoo_c_2D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 2);
D.fn = @(xh) branin_hoo_orig(xh);
PB.f_star      = 0.39788735772973816;
PB.f_star_justification = 'A';   % analytic reduction, 1-D root to machine precision
PB.n_minimizers_documented = 3;
PB.x0_idx      = [5.7666015625 3.5459533607681757];   % protocol starting point: Halton point 1005 (index representation)
PB.x0_alt_idx  = [-1.94580078125 3.9780521262002742;7.54638671875 4.9039780521262006;-0.01220703125 3.8545953360768181;3.85498046875 4.1220850480109741;-3.79150390625 4.554183813443073;8.51318359375 3.5048010973936896;0.95458984375 4.4307270233196165;4.82177734375 4.6433470507544579;-2.86865234375 3.5939643347050754];   % 9 alternative starting points (Halton indices 1005 + 63 k)
PB.core_pidx   = 5;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('branin_hoo_c_2D: n = 2, pattern C, m = 0, f* = 0.3978873577 (A)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 2 || ~all(isfinite(z)), error('branin_hoo_c_2D: z must be a finite numeric vector with 2 entries.'); end
xh = zeros(2, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('branin_hoo_c_2D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('branin_hoo_c_2D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
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
if ~iscell(xr) || numel(xr) ~= 2, error('branin_hoo_c_2D: xr must be a cell array with 2 entries.'); end
z = zeros(1, 2);
for i = 1:2
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('branin_hoo_c_2D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('branin_hoo_c_2D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('branin_hoo_c_2D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('branin_hoo_c_2D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = branin_hoo_orig(x)
%
% Purpose:
%
%    Function brannin_hoo is the function described in
%    Garcia-Palomares(2012).
%
%    dim = 2
%    Minimum global value =  5/(4*pi)
%    Global minima = ([-pi 12.275]', [pi 2.275]', [3*pi 2.475]')
%    Local minima (equal to global minima)
%    Search domain: -5 <= x <= 20 (Garcia-Palomares(2012))
%
%                    -5 <= x(1) <= 10 
%                     0 <= x(2) <= 15 (Huyer and Neumaier (1999))
%                                     (Huyer and Neumaier (2008))
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
f = (x(2)-1.275*(x(1)/pi)^2+5*x(1)/pi-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end
%
% End of branin_hoo.

