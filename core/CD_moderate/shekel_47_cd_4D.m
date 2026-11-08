function PB = shekel_47_cd_4D()
%SHEKEL_47_CD_4D  DFO Benchmark Suite v2.1.0 core: shekel_47, dimension 4, pattern CD, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 151 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = -10.40294056681866 (L: literature value at a refined point); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:D x3:C x4:D
%
lb_orig = [1.3487915039062501 1.5352893359752071 1.8390336 1.9793835901707624];
ub_orig = [10 10 10 10];
lb_C_work = [1.3487915039062504 1.5352893359752071 183.90336000000002 1.9793835901707624];
ub_C_work = [10 10 1000 10];
T_break = [NaN 4.0182711307558128 NaN 4.0112730806608354];   % breakpoints of T (NaN for continuous)
T_image = [NaN 4.0006893663807626 NaN 3.9996061601069748];   % their images
m = 151;
ProblemData = cell(1, 4);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = {linspace(lb_orig(2), ub_orig(2), m)'};
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = {linspace(lb_orig(4), ub_orig(4), m)'};
PB.name        = 'shekel_47_cd_4D';
PB.base_function = 'shekel_47';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CD adaptation of shekel_47 (F = f o T).';
PB.nvar        = 4;
PB.nobj        = 1;
PB.var_types   = {'C','D','C','D'};
PB.type_counts = struct('C', 2, 'D', 2, 'K', 0);
PB.proportion  = [2 2 0];
PB.pattern     = 'CD';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100];
PB.cardinality = [Inf 151 Inf 151];
PB.n_finite    = 2;
PB.finite_part_size = 22801;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'shekel_47_cd_4D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 4);
D.cw = [4.0005729162446064 NaN 399.94897114568784 NaN];   % anchors of the continuous maps (NaN for finite)
D.cx = [4.0005729162446064 NaN 3.9994897114568779 NaN];
D.fn = @(xh) shekel_47_orig(xh);
PB.f_star      = -10.402940566818662;
PB.f_star_justification = 'L';   % literature value at a refined point
PB.n_minimizers_documented = 1;
PB.x0_idx      = [2.1640665233135228 83 814.84399431680015 114];   % protocol starting point: Halton point 1048 (index representation)
PB.x0_alt_idx  = [9.2607609927654266 72 442.70392647680006 121;4.9013629615306851 77 979.36907694080003 128;7.1317526519298555 79 607.2290091008 113;2.7216639459133152 69 197.48320808960003 120;8.1962568223476406 78 765.48646952960007 126;3.8368587911129 82 362.00829071359999 109;6.0672484815120695 70 930.01155215360006 115;1.6318144381046298 74 526.53337333759998 122;9.5395597040653222 83 311.08386037759999 129];   % 9 alternative starting points (Halton indices 1048 + 63 k)
PB.core_pidx   = 48;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('shekel_47_cd_4D: n = 4, pattern CD, m = 151, f* = -10.40294057 (L)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 4 || ~all(isfinite(z)), error('shekel_47_cd_4D: z must be a finite numeric vector with 4 entries.'); end
xh = zeros(4, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('shekel_47_cd_4D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('shekel_47_cd_4D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
k = check_idx(z(4), D.m, D.name, 4); g = D.lb(4) + (k - 1) * (D.ub(4) - D.lb(4)) / (D.m - 1);
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
if ~iscell(xr) || numel(xr) ~= 4, error('shekel_47_cd_4D: xr must be a cell array with 4 entries.'); end
z = zeros(1, 4);
for i = 1:4
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_47_cd_4D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('shekel_47_cd_4D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('shekel_47_cd_4D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('shekel_47_cd_4D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = shekel_47_orig(x)
%
% Purpose:
%
%    Function shekel_4 is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value  
%    Global minima = 4*ones(4,1)
%    Local minima (m local minima at x = A(:,1:4))
%    Search domain: 0 <= x <= 10 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
%    Cases considered: m = 5 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 7 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
%                      m = 10 Brachetti et al. (1997)
%                            Huyer and Neumaier (1999)
%                            Huyer and Neumaier (2008)
%                            Jones et al. (1993)
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
m = 7;
A = [4 4 4 4 0.1
    1 1 1 1 0.2
    8 8 8 8 0.2
    6 6 6 6 0.4
    3 7 3 7 0.4
    2 9 2 9 0.6
    5 5 3 3 0.3
    8 1 8 1 0.7
    6 2 6 2 0.5
    7 3.6 7 3.6 0.5];
f = -sum(1./(sum((repmat(x',m,1)-A(1:m,1:4)).^2,2)+A(1:m,5)),1);
end
%
% End of shekel_4.

