function PB = transistor_c_9D()
%TRANSISTOR_C_9D  DFO Benchmark Suite v2 core: transistor, dimension 9, pattern C, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the original interval, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C
%
lb_orig = [-10 -10 -10 -10 -10 -10 -10 -10 -10];
ub_orig = [10 10 10 10 10 10 10 10 10];
lb_C_work = [-10 -1000 -0.10000000000000001 -10 -1000 -0.10000000000000001 -10 -1000 -0.10000000000000001];
ub_C_work = [10 1000 0.10000000000000001 10 1000 0.10000000000000001 10 1000 0.10000000000000001];
T_break = [NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 9);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = [lb_C_work(6); ub_C_work(6)];
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = [lb_C_work(8); ub_C_work(8)];
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
PB.name        = 'transistor_c_9D';
PB.base_function = 'transistor';
PB.suite_version = 2;
PB.core = true;
PB.description = 'DFO Benchmark Suite v2 core: C adaptation of transistor (F = f o T).';
PB.nvar        = 9;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 9, 'D', 0, 'K', 0);
PB.proportion  = [9 0 0];
PB.pattern     = 'C';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01 1 100 0.01 1 100 0.01];
PB.cardinality = [Inf Inf Inf Inf Inf Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'transistor_c_9D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 9);
D.fn = @(xh) transistor_orig(xh);
PB.f_star      = 1.6004932665473596e-24;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-2.177734375 -914.0374942844079 0.00022400000000000198 5.943356934610577 235.16153268219387 0.042558033682294022 0.18522287807856763 884.23968508528924 -0.065184515492726219];   % protocol starting point: Halton point 1062 (index representation)
PB.x0_alt_idx  = [2.978515625 -790.58070416095109 -0.098656000000000008 6.8179925031236976 -391.43501126972205 0.017705962676376891 -5.4203134541013638 -566.99227292608259 0.083397715131092287;-7.099609375 -930.49839963420209 0.040543999999999997 4.8937942523948337 -837.71600300525915 -0.022439690487027764 7.7325463057195201 81.207173057297041 0.032358017588559218;5.556640625 -894.83310470964796 -0.059872000000000002 5.7684298209079525 519.15852742299035 -0.04729176149294493 2.1270099735395913 729.40661904067656 -0.018681679953973879;-4.638671875 -837.21993598536812 0.07932800000000001 6.6430653894210749 72.877535687453246 -0.072143832498862082 -3.4785263586403428 -616.85376877095791 -0.069721377496506948;0.517578125 -977.13763145861913 -0.019552 4.3190337359433535 -570.11133119322449 0.086527082385070564 -9.0840626908202733 31.345677212421606 0.078860853127311586;-9.560546875 -853.68084133516231 -0.072352 5.1936693044564741 967.0787514514036 0.061675011379153405 5.2411968247506593 679.54512319580112 0.02782115558477849;9.658203125 -825.33150434385004 0.020768000000000009 6.0683048729695912 340.48220749948769 0.036822940373236246 -1.5367392631793209 -771.68683481557082 -0.023218541957754579;-0.595703125 -965.24919981710104 -0.032032000000000005 6.9429404414827118 -105.79878423604953 -0.0033227127901684084 -7.1422755953592514 -123.4873888321913 -0.074258239500287662;4.560546875 -907.63603109282121 0.061088000000000003 5.0187421907538514 -732.39532818796533 -0.028174783796085567 7.182983920211683 524.71205715118799 0.065644776855428633];   % 9 alternative starting points (Halton indices 1062 + 63 k)
PB.core_pidx   = 62;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('transistor_c_9D: n = 9, pattern C, m = 0, f* = 1.600493267e-24 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 9 || ~all(isfinite(z)), error('transistor_c_9D: z must be a finite numeric vector with 9 entries.'); end
xh = zeros(9, 1);
v = z(1); if v < D.lbw(1) || v > D.ubw(1), error('transistor_c_9D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1));
v = z(2); if v < D.lbw(2) || v > D.ubw(2), error('transistor_c_9D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2));
v = z(3); if v < D.lbw(3) || v > D.ubw(3), error('transistor_c_9D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3));
v = z(4); if v < D.lbw(4) || v > D.ubw(4), error('transistor_c_9D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4));
v = z(5); if v < D.lbw(5) || v > D.ubw(5), error('transistor_c_9D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5));
v = z(6); if v < D.lbw(6) || v > D.ubw(6), error('transistor_c_9D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6));
v = z(7); if v < D.lbw(7) || v > D.ubw(7), error('transistor_c_9D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7));
v = z(8); if v < D.lbw(8) || v > D.ubw(8), error('transistor_c_9D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8));
v = z(9); if v < D.lbw(9) || v > D.ubw(9), error('transistor_c_9D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9));
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
xr = cell(1, 9);
for i = 1:9
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 9, error('transistor_c_9D: xr must be a cell array with 9 entries.'); end
z = zeros(1, 9);
for i = 1:9
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('transistor_c_9D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v), error('transistor_c_9D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('transistor_c_9D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('transistor_c_9D: x%d unknown label %s.', i, lab); end
            z(i) = k;
    end
end
end

% -------------------------------------------------------------------------
% Embedded original continuous function (verbatim; only the top-level name is renamed)
% -------------------------------------------------------------------------
function [f] = transistor_orig(x)
%
% Purpose:
%
%    Function transistor is the function described in
%    Montaz et al. (2005).
%
%    dim = 9
%    Minimum global value = 0
%    Global minima near [0.9 0.45 1 2 8 8 5 1 2]'
%    Local minima
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
g = [0.485 0.752 0.869 0.982
    0.369 1.254 0.703 1.455
    5.2095 10.0677 22.9274 20.2153
    23.3037 101.779 111.461 191.267
    28.5132 111.8467 134.3884 211.4823]';
alfa = (1-x(1)*x(2))*x(3).*(exp(x(5).*(g(:,1)-g(:,3).*x(7)*0.001-...
       g(:,5).*x(8)*0.001))-1)-g(:,5)+g(:,4).*x(2);
beta = (1-x(1)*x(2))*x(4).*(exp(x(6).*(g(:,1)-g(:,2)-g(:,3).*x(7)*0.001+...
       g(:,4).*x(9)*0.001))-1)-g(:,5).*x(1)+g(:,4);
f = (x(1)*x(3)-x(2)*x(4))^2 + sum(alfa.^2+beta.^2);
end
%
% End of transistor.

