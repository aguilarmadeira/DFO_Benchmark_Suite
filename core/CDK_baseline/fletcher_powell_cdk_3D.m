function PB = fletcher_powell_cdk_3D()
%FLETCHER_POWELL_CDK_3D  DFO Benchmark Suite v2.1.0 core: fletcher_powell, dimension 3, pattern CDK, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:K x3:D
%
lb_orig = [-10 -10 -10];
ub_orig = [7.1000732421875004 7.351917898694305 7.1731327999999994];
lb_C_work = [-10 -10 -10];
ub_C_work = [7.1000732421875004 7.351917898694305 7.1731327999999994];
T_break = [NaN -0.051567071415265175 -0.039582976000001935];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0 0];   % their images
m = 151;
ProblemData = cell(1, 3);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = arrayfun(@(j) sprintf('k2_%d', j), (1:m)', 'UniformOutput', false);
ProblemData{3} = {linspace(lb_orig(3), ub_orig(3), m)'};
PB.name        = 'fletcher_powell_cdk_3D';
PB.base_function = 'fletcher_powell';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CDK adaptation of fletcher_powell (F = f o T).';
PB.nvar        = 3;
PB.nobj        = 1;
PB.var_types   = {'C','K','D'};
PB.type_counts = struct('C', 1, 'D', 1, 'K', 1);
PB.proportion  = [1 1 1];
PB.pattern     = 'CDK';
PB.strategy    = 'baseline';
PB.scale_factors_C = 1;
PB.cardinality = [Inf 151 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 151;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'fletcher_powell_cdk_3D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 3);
D.cw = [1 NaN NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [1 NaN NaN];
g = linspace(lb_orig(2), ub_orig(2), m); D.E{2} = g([36 71 104 99 74 107 147 131 93 5 144 41 51 45 16 135 43 129 112 150 35 30 122 138 21 102 62 125 97 67 53 19 87 26 79 80 92 47 143 64 54 109 126 60 39 4 10 142 146 25 57 96 130 134 103 66 91 90 120 94 95 85 7 15 33 48 58 117 118 38 44 2 27 115 68 111 29 110 77 32 59 139 6 56 3 141 22 34 61 69 55 72 116 65 84 151 124 73 88 14 89 42 70 23 106 105 37 8 136 52 121 28 145 82 40 13 148 11 98 137 127 149 108 132 75 119 100 31 114 17 49 9 123 20 1 50 18 128 113 24 83 12 101 63 78 81 140 86 76 133 46])';
D.fn = @(xh) fletcher_powell_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-5.7416809797286987 31 26];   % protocol starting point: Halton point 1020 (index representation)
PB.x0_alt_idx  = [4.7037250876426686 20 101;-3.712717211246491 25 56;0.49550393819808924 34 132;-7.9209383606910713 22 85;5.8058782458305345 26 17;-2.6105640530586252 30 115;1.5976570963859551 18 47;-6.8187852025032045 27 146;3.6516698002815229 31 77];   % 9 alternative starting points (Halton indices 1020 + 63 k)
PB.core_pidx   = 20;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('fletcher_powell_cdk_3D: n = 3, pattern CDK, m = 151, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 3 || ~all(isfinite(z)), error('fletcher_powell_cdk_3D: z must be a finite numeric vector with 3 entries.'); end
xh = zeros(3, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('fletcher_powell_cdk_3D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.E{2}(k);
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
if ~iscell(xr) || numel(xr) ~= 3, error('fletcher_powell_cdk_3D: xr must be a cell array with 3 entries.'); end
z = zeros(1, 3);
for i = 1:3
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('fletcher_powell_cdk_3D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('fletcher_powell_cdk_3D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('fletcher_powell_cdk_3D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('fletcher_powell_cdk_3D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = fletcher_powell_orig(x)
%
% Purpose:
%
%    Function fletcher_powell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 3
%    Minimum global value = 0
%    Global minima = [1 0 0]'
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
if x(1) == 0
    f = +inf;
else
   if x(1)>0
      teta = 1/(2*pi)*atan(x(2)/x(1));
   else
      teta = 1/(2*pi)*atan(x(2)/x(1))+0.5;
   end
   f = 100* ((x(3)-10*teta)^2+(norm(x,2)-1)^2)+x(3)^2;
end
end
%
% End of fletcher_powell.

