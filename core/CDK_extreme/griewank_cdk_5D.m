function PB = griewank_cdk_5D()
%GRIEWANK_CDK_5D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 5, pattern CDK, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 41 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D x4:C x5:D
%
lb_orig = [-600 -600 -379.98796800000002 -600 -600];
ub_orig = [363.00439453125 431.78174058832496 600 415.61682632236568 405.60016392322927];
lb_C_work = [-600 -600 -379.98796800000002 -60000000000 -600];
ub_C_work = [363.00439453125 431.78174058832496 600 41561682632.236572 405.60016392322927];
T_break = [NaN -6.7254991617130599 12.007219200000009 NaN 3.3600983539375875];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 0 NaN 0];   % their images
m = 41;
ProblemData = cell(1, 5);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = {linspace(lb_orig(5), ub_orig(5), m)'};
PB.name        = 'griewank_cdk_5D';
PB.base_function = 'griewank';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CDK adaptation of griewank (F = f o T).';
PB.nvar        = 5;
PB.nobj        = 1;
PB.var_types   = {'C','K','D','C','D'};
PB.type_counts = struct('C', 2, 'D', 2, 'K', 1);
PB.proportion  = [2 2 1];
PB.pattern     = 'CDK';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 100000000];
PB.cardinality = [Inf 41 41 Inf 41];
PB.n_finite    = 3;
PB.finite_part_size = 68921;
PB.M_K = 41;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'griewank_cdk_5D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 5);
D.cw = [0 NaN NaN 0 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [0 NaN NaN 0 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([15 16 18 33 40 1 32 11 29 22 30 3 19 12 25 41 14 10 21 37 8 4 7 34 39 35 9 20 27 2 31 5 13 26 23 36 17 6 24 38 28])';
D.fn = @(xh) griewank_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [362.06396055221558 13 32 -31193875105.142406 3];   % protocol starting point: Halton point 1023 (index representation)
PB.x0_alt_idx  = [-133.07452940940857 11 13 -42995503366.031197 34;115.20004105567932 12 40 -38554030364.621429 21;-370.06389212608337 14 21 -34112557363.21167 9;239.33732628822327 11 1 -43883797966.313148 39;-251.56921076774597 12 29 -39442324964.903381 26;-3.2946403026580811 13 9 -37031239621.280945 18;-488.55857348442078 10 37 -32589766619.871185 5;301.40596890449524 12 17 -42361007222.972656 36;-192.32187008857727 14 7 -37919534221.562897 23];   % 9 alternative starting points (Halton indices 1023 + 63 k)
PB.core_pidx   = 23;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('griewank_cdk_5D: n = 5, pattern CDK, m = 41, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 5 || ~all(isfinite(z)), error('griewank_cdk_5D: z must be a finite numeric vector with 5 entries.'); end
xh = zeros(5, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('griewank_cdk_5D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
k = check_idx(z(3), D.m, D.name, 3); g = D.lb(3) + (k - 1) * (D.ub(3) - D.lb(3)) / (D.m - 1);
xh(3) = T_map(g, D.lb(3), D.ub(3), D.Tb(3), D.Ti(3));
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('griewank_cdk_5D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
k = check_idx(z(5), D.m, D.name, 5); g = D.lb(5) + (k - 1) * (D.ub(5) - D.lb(5)) / (D.m - 1);
xh(5) = T_map(g, D.lb(5), D.ub(5), D.Tb(5), D.Ti(5));
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
if ~iscell(xr) || numel(xr) ~= 5, error('griewank_cdk_5D: xr must be a cell array with 5 entries.'); end
z = zeros(1, 5);
for i = 1:5
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_cdk_5D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_cdk_5D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('griewank_cdk_5D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('griewank_cdk_5D: x%d unknown label %s.', i, lab); end
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

