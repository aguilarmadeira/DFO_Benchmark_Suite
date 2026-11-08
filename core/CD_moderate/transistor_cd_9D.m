function PB = transistor_cd_9D()
%TRANSISTOR_CD_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern CD, strategy moderate.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 25 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:D x3:C x4:D x5:C x6:D x7:C x8:D x9:C
%
lb_orig = [-10 -10 -10 -10 -7.3500853766819212 -6.9926613213823039 -10 -10 -10];
ub_orig = [6.8211669921875 7.0331250317532898 6.7027327999999997 6.1800083298625577 10 10 6.7876079069934505 6.4033133570184386 7.3164625626695159];
lb_C_work = [-10 -10 -1000 -10 -0.07350085376681921 -6.9926613213823039 -10 -10 -1000];
ub_C_work = [6.8211669921875 7.0331250317532898 670.27327999999989 6.1800083298625577 0.099999999999999992 10 6.7876079069934505 6.4033133570184386 731.64625626695147];
T_break = [NaN 0.64570314484580571 NaN 2.1350062473969214 NaN 7.8759173348272125 NaN 0.93554223801229419 NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN 0.44998747198158512 NaN 2.000068541623683 NaN 7.9996926842185427 NaN 0.99998772345678977 NaN];   % their images
m = 25;
ProblemData = cell(1, 9);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = {linspace(lb_orig(2), ub_orig(2), m)'};
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = {linspace(lb_orig(4), ub_orig(4), m)'};
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = {linspace(lb_orig(6), ub_orig(6), m)'};
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = {linspace(lb_orig(8), ub_orig(8), m)'};
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
PB.name        = 'transistor_cd_9D';
PB.base_function = 'transistor';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: CD adaptation of transistor (F = f o T).';
PB.nvar        = 9;
PB.nobj        = 1;
PB.var_types   = {'C','D','C','D','C','D','C','D','C'};
PB.type_counts = struct('C', 5, 'D', 4, 'K', 0);
PB.proportion  = [5 4 0];
PB.pattern     = 'CD';
PB.strategy    = 'moderate';
PB.scale_factors_C = [1 100 0.01 1 100];
PB.cardinality = [Inf 25 Inf 25 Inf 25 Inf 25 Inf];
PB.n_finite    = 4;
PB.finite_part_size = 390625;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'transistor_cd_9D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 9);
D.cw = [0.89999995261685584 NaN 100.00064824653873 NaN 0.079999714405078512 NaN 5.0000312759296719 NaN 200.00524834861903];   % anchors of the continuous maps (NaN for finite)
D.cx = [0.8999999526168565 NaN 1.0000064824653869 NaN 7.9999714405078528 NaN 5.0000312759296719 NaN 2.0000524834861899];
D.fn = @(xh) transistor_orig(xh);
PB.f_star      = 1.6004932665473596e-24;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-3.4210181832313538 2 -162.99265392640007 20 0.03364993646332759 18 -1.4507235938738816 24 -698.55948296429665];   % protocol starting point: Halton point 1062 (index representation)
PB.x0_alt_idx  = [0.91568893194198608 4 -988.77576355839994 21 -0.020707581208170241 15 -6.1559008965260258 6 587.89983407334398;-7.5606022477149963 2 173.73443932159989 19 -0.059422647751179489 10 4.884351728651227 14 145.98632822071932;3.084042489528656 2 -664.87636910080005 20 0.058286796990697123 7 0.17917442599908284 22 -295.92717763190535;-5.490810215473175 3 497.6338337791999 21 0.019571730447687868 4 -4.5260028766530622 6 -737.8406834845299;-1.1541031002998352 1 -328.14927585280009 18 -0.036207828235490742 23 -9.2311801793052055 13 548.61863355311061;-9.6303942799568176 3 -769.10142177279999 19 0.097144067634874418 20 2.7931618163613834 21 106.70512770048595;6.5336958765983582 3 8.5778173951999861 20 0.042786549963376594 17 -2.8961048567800978 4 -335.20837815213849;-2.0904375910758972 1 -432.37432852480003 21 0.0040714834203673528 13 -7.6012821594322411 12 -777.12188400476316;2.2462695240974426 2 345.30491064319995 19 -0.050286034251130471 10 4.4230598362343478 19 434.19078855938778];   % 9 alternative starting points (Halton indices 1062 + 63 k)
PB.core_pidx   = 62;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('transistor_cd_9D: n = 9, pattern CD, m = 25, f* = 1.600493267e-24 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 9 || ~all(isfinite(z)), error('transistor_cd_9D: z must be a finite numeric vector with 9 entries.'); end
xh = zeros(9, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('transistor_cd_9D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
k = check_idx(z(2), D.m, D.name, 2); g = D.lb(2) + (k - 1) * (D.ub(2) - D.lb(2)) / (D.m - 1);
xh(2) = T_map(g, D.lb(2), D.ub(2), D.Tb(2), D.Ti(2));
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('transistor_cd_9D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
k = check_idx(z(4), D.m, D.name, 4); g = D.lb(4) + (k - 1) * (D.ub(4) - D.lb(4)) / (D.m - 1);
xh(4) = T_map(g, D.lb(4), D.ub(4), D.Tb(4), D.Ti(4));
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('transistor_cd_9D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
k = check_idx(z(6), D.m, D.name, 6); g = D.lb(6) + (k - 1) * (D.ub(6) - D.lb(6)) / (D.m - 1);
xh(6) = T_map(g, D.lb(6), D.ub(6), D.Tb(6), D.Ti(6));
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('transistor_cd_9D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
k = check_idx(z(8), D.m, D.name, 8); g = D.lb(8) + (k - 1) * (D.ub(8) - D.lb(8)) / (D.m - 1);
xh(8) = T_map(g, D.lb(8), D.ub(8), D.Tb(8), D.Ti(8));
v = z(9); tl = 1e-12 * (D.ubw(9) - D.lbw(9)); if v < D.lbw(9) - tl || v > D.ubw(9) + tl, error('transistor_cd_9D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
v = min(max(v, D.lbw(9)), D.ubw(9));
if abs(v - D.cw(9)) <= 1e-12 * (D.ubw(9) - D.lbw(9)), xh(9) = D.cx(9); else, xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9)); end
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
if ~iscell(xr) || numel(xr) ~= 9, error('transistor_cd_9D: xr must be a cell array with 9 entries.'); end
z = zeros(1, 9);
for i = 1:9
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('transistor_cd_9D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('transistor_cd_9D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('transistor_cd_9D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('transistor_cd_9D: x%d unknown label %s.', i, lab); end
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

