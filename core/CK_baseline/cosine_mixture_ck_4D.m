function PB = cosine_mixture_ck_4D()
%COSINE_MIXTURE_CK_4D  DFO Benchmark Suite v2.1.0 core: cosine_mixture, dimension 4, pattern CK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -4.4 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:C x4:K
%
lb_orig = [-0.66625732421874995 -0.65395722196819595 -0.62211327999999999 -0.69840899625156183];
ub_orig = [1 1 1 1];
lb_C_work = [-0.66625732421874995 -0.65395722196819595 -0.62211327999999999 -0.69840899625156183];
ub_C_work = [1 1 1 1];
T_break = [NaN 1 NaN 1];   % breakpoints of T (NaN for continuous)
T_image = [NaN 1 NaN 1];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = arrayfun(@(j) sprintf('k4_%d', j), (1:m)', 'UniformOutput', false);
PB.name        = 'cosine_mixture_ck_4D';
PB.base_function = 'cosine_mixture';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CK adaptation of cosine_mixture (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','K','C','K'};
PB.type_counts = struct('C', 2, 'D', 0, 'K', 2);
PB.proportion  = [2 0 2];
PB.pattern     = 'CK';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1];
PB.cardinality = [Inf 151 Inf 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 22801;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'cosine_mixture_ck_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [1 NaN 1 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [1 NaN 1 NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([38 121 109 75 118 138 137 65 148 147 56 44 35 108 62 26 77 93 17 60 55 33 3 94 36 90 53 98 66 23 103 101 115 119 6 8 124 145 122 149 82 80 37 68 2 150 16 54 141 72 96 81 15 89 128 144 88 70 7 19 39 105 50 58 112 21 29 84 135 120 79 92 86 146 12 125 43 110 116 113 57 74 130 42 69 52 73 85 18 123 102 97 27 71 134 142 91 48 139 106 31 76 127 46 107 117 131 59 95 63 104 9 34 20 64 129 67 49 4 1 51 22 11 13 30 140 100 136 151 83 99 5 132 87 126 32 61 14 111 47 24 40 45 41 114 78 143 28 25 10 133])';
g = linspace(lb_orig(4), ub_orig(4), m); D.E{4} = g([145 28 117 11 108 103 120 52 29 74 101 112 144 99 82 26 37 15 128 151 127 134 122 92 69 106 138 19 73 115 30 50 102 91 131 148 14 60 66 27 31 123 141 44 70 65 90 54 89 32 42 98 17 20 100 133 81 136 71 116 135 113 53 61 12 147 33 139 51 109 9 93 104 67 16 114 96 118 10 129 62 87 105 47 126 58 94 97 25 18 83 142 34 79 95 48 3 137 55 76 40 86 111 56 5 78 4 121 107 77 7 140 49 35 149 22 84 150 23 21 72 132 36 38 57 75 1 63 143 119 8 45 43 2 41 88 125 85 46 130 124 6 39 64 13 24 146 59 80 68 110])';
D.fn = @(xh) cosine_mixture_orig(xh);
PB.f_star      = -4.4000000000000004;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [0.26938521623611456 59 0.74876709519360007 37];   % protocol starting point: Halton point 1009 (index representation)
PB.x0_alt_idx  = [-0.58733790993690482 63 0.32052918927360019 41;0.93572542548179616 53 -0.4814436163584 27;0.096087945699691812 62 0.64754722652160024 33;0.52566991209983827 64 -0.16688340910080002 40;-0.32373079419136042 54 0.96210743377920016 25;0.73069766879081732 58 0.16013462814720003 29;-0.10893981099128724 68 -0.57954902753280002 35;0.32064215540885932 55 0.48715266539520008 42;-0.53364016413688653 59 -0.25253099028479997 27];   % 9 alternative starting points (Halton indices 1009 + 63 k)
PB.core_pidx   = 9;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('cosine_mixture_ck_4D: n = 4, pattern CK, m = 151, f* = -4.4 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('cosine_mixture_ck_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('cosine_mixture_ck_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('cosine_mixture_ck_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
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
if ~iscell(xr) || numel(xr) ~= 4, error('cosine_mixture_ck_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('cosine_mixture_ck_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('cosine_mixture_ck_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('cosine_mixture_ck_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('cosine_mixture_ck_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = cosine_mixture_orig(x)
%
% Purpose:
%
%    Function cosine_mixture is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value: -2 if n = 2
%                          -4 if n = 4
%    Global minima
%    Local minima
%    Search domain: -1 <= x <= 1 (Brachetti et al. (1997))
%    Cases considered: n = 2 Brachetti et al. (1997)
%                      n = 4 Brachetti et al. (1997)
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
f = 0.1*sum(cos(5*pi.*x),1) - sum(x.^2,1); 
end
%
% End of cosine_mixture.

