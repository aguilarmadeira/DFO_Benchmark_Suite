function PB = miele_cantrell_c_4D()
%MIELE_CANTRELL_C_4D  DFO Benchmark Suite v2 core: miele_cantrell, dimension 4, pattern C, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C
%
lb_orig = [-10 -10 -10 -10];
ub_orig = [10 10 10 10];
lb_C_work = [-10 -10 -1000000000 -1000000000];
ub_C_work = [10 10 1000000000 1000000000];
T_break = [NaN NaN NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
PB.name        = 'miele_cantrell_c_4D';
PB.base_function = 'miele_cantrell';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of miele_cantrell (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C'};
PB.type_counts = struct('C', 4, 'D', 0, 'K', 0);
PB.proportion  = [4 0 0];
PB.pattern     = 'C';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 1 100000000 100000000];
PB.cardinality = [Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'miele_cantrell_c_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.fn = @(xh) miele_cantrell_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-8.740234375 -5.0251486053955192 -93760000 -140358184.0899626];   % protocol starting point: Halton point 1032 (index representation)
PB.x0_alt_idx  = [7.666015625 -4.4490169181527204 -621760000 -52894627.23865068;-2.412109375 -3.8728852309099215 309440000.00000024 34568929.612661481;2.744140625 -5.4915409236396897 -218560000 122032486.46397328;-7.451171875 -4.2569730224051208 697280000.00000024 -70387338.608912945;5.205078125 -3.6808413351623228 169280000 17076218.242399096;-4.873046875 -5.0800182898948334 -819520000 64556434.818825245;0.283203125 -4.7965249199817102 572480000.00000024 -127863390.25406086;-9.970703125 -3.561957018747143 -416320000 -40399833.402748942;9.716796875 -4.9611339734796536 975680000 47063723.44856286];   % 9 alternative starting points (Halton indices 1032 + 63 k)
PB.core_pidx   = 32;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('miele_cantrell_c_4D: n = 4, pattern C, m = 0, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('miele_cantrell_c_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('miele_cantrell_c_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('miele_cantrell_c_4D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('miele_cantrell_c_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('miele_cantrell_c_4D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
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
if ~iscell(xr) || numel(xr) ~= 4, error('miele_cantrell_c_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('miele_cantrell_c_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('miele_cantrell_c_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('miele_cantrell_c_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('miele_cantrell_c_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = miele_cantrell_orig(x)
%
% Purpose:
%
%    Function miele_cantrell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [0 1 1 1]'
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
f = (exp(x(1))-x(2))^4 + 100*(x(2)-x(3))^6 + tan(x(3)-x(4))^4 + x(1)^8;
end
%
% End of miele_cantrell.

