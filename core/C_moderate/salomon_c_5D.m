function PB = salomon_c_5D()
%SALOMON_C_5D  DFO Benchmark Suite v2 core: salomon, dimension 5, pattern C, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C
%
lb_orig = [-100 -100 -100 -100 -100];
ub_orig = [100 100 100 100 100];
lb_C_work = [-100 -10000 -1 -100 -10000];
ub_C_work = [100 10000 1 100 10000];
T_break = [NaN NaN NaN NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN NaN NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 5);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
PB.name        = 'salomon_c_5D';
PB.base_function = 'salomon';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of salomon (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C'};
PB.type_counts = struct('C', 5, 'D', 0, 'K', 0);
PB.proportion  = [5 0 0];
PB.pattern     = 'C';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01 1 100];
PB.cardinality = [Inf Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'salomon_c_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
D.fn = @(xh) salomon_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-43.65234375 2382.2588020118874 0.066240000000000077 75.760099958350679 5657.400450788884];   % protocol starting point: Halton point 1042 (index representation)
PB.x0_alt_idx  = [7.91015625 2958.3904892546852 -0.84575999999999996 84.506455643481871 -608.56498873027886;-92.87109375 1559.2135345221777 0.46944000000000008 93.252811328613063 -5071.3749060856499;89.94140625 1915.8664837677188 -0.44255999999999995 74.010828821324424 8497.3703981968465;-12.01171875 3150.4343850022869 0.85728000000000004 82.757184506455673 4034.5604808414737;39.55078125 1751.2574302697776 -0.054720000000000102 91.503540191586836 -2395.3281879653032;-61.23046875 2327.3891175125755 -0.65952 96.251561849229461 -6858.1381053206751;65.33203125 2610.8824874256989 0.34848000000000035 77.009579341940849 6710.6071989618213;-37.20703125 1211.7055326931877 -0.2563200000000001 85.755935027072013 2247.7972816064466;14.35546875 2446.273433927754 0.75168000000000013 94.502290712203205 -4018.1681579127116];   % 9 alternative starting points (Halton indices 1042 + 63 k)
PB.core_pidx   = 42;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('salomon_c_5D: n = 5, pattern C, m = 0, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('salomon_c_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('salomon_c_5D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('salomon_c_5D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('salomon_c_5D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('salomon_c_5D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('salomon_c_5D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
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
if ~iscell(xr) || numel(xr) ~= 5, error('salomon_c_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('salomon_c_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('salomon_c_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('salomon_c_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('salomon_c_5D: x%d unknown label %s.', i, lab); end
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

