function PB = griewank_c_10D()
%GRIEWANK_C_10D  DFO Benchmark Suite v2.1.0 core: griewank, dimension 10, pattern C, strategy baseline.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 0 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C x10:C
%
lb_orig = [-270.0029296875 -256.74338261443881 -400 -400 -400 -400 -400 -400 -278.2615270814498 -287.59161917257779];
ub_orig = [400 400 264.52531199999999 285.07788421491045 275.49101837306193 285.00822800322118 291.6565654146861 262.337520430322 400 400];
lb_C_work = [-270.0029296875 -256.74338261443881 -400 -400 -400 -400 -400 -400 -278.2615270814498 -287.59161917257779];
ub_C_work = [400 400 264.52531199999999 285.07788421491045 275.49101837306193 285.00822800322118 291.6565654146861 262.337520430322 400 400];
T_break = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % breakpoints of T (NaN for continuous)
T_image = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];   % their images
m = 0;
ProblemData = cell(1, 10);
ProblemData{1} = [lb_C_work(1); ub_C_work(1)];
ProblemData{2} = [lb_C_work(2); ub_C_work(2)];
ProblemData{3} = [lb_C_work(3); ub_C_work(3)];
ProblemData{4} = [lb_C_work(4); ub_C_work(4)];
ProblemData{5} = [lb_C_work(5); ub_C_work(5)];
ProblemData{6} = [lb_C_work(6); ub_C_work(6)];
ProblemData{7} = [lb_C_work(7); ub_C_work(7)];
ProblemData{8} = [lb_C_work(8); ub_C_work(8)];
ProblemData{9} = [lb_C_work(9); ub_C_work(9)];
ProblemData{10} = [lb_C_work(10); ub_C_work(10)];
PB.name        = 'griewank_c_10D';
PB.base_function = 'griewank';
PB.suite_version = 2;
PB.core = true;
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: C adaptation of griewank (F = f o T).';
PB.nvar        = 10;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 10, 'D', 0, 'K', 0);
PB.proportion  = [10 0 0];
PB.pattern     = 'C';
PB.strategy    = 'baseline';
PB.scale_factors_C = [1 1 1 1 1 1 1 1 1 1];
PB.cardinality = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];
PB.n_finite    = 0;
PB.finite_part_size = 1;
PB.M_K = 1;
PB.domain_size = Inf;
PB.lb_orig = lb_orig(:); PB.ub_orig = ub_orig(:);
PB.lb_C_work = lb_C_work(:); PB.ub_C_work = ub_C_work(:);
PB.m = m;
PB.ProblemData = ProblemData;
D = struct('name', 'griewank_c_10D', 'm', m, 'lb', lb_orig(:), 'ub', ub_orig(:), 'lbw', lb_C_work(:), 'ubw', ub_C_work(:), 'Tb', T_break(:), 'Ti', T_image(:));
D.E = cell(1, 10);
D.cw = [0 0 5.6843418860808015e-14 -5.6843418860808015e-14 0 -5.6843418860808015e-14 5.6843418860808015e-14 5.6843418860808015e-14 0 0];   % anchors of the continuous maps (NaN for finite)
D.cx = [0 0 0 0 0 0 0 0 0 0];
D.fn = @(xh) griewank_orig(xh);
PB.f_star      = 0;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [64.344235420227051 312.01380379239561 -24.463455682559982 -303.55838198057484 240.47307677445849 23.414280213188135 -296.66681884502754 150.61212151825282 43.615850856274449 -116.68750685045856];   % protocol starting point: Halton point 1022 (index representation)
PB.x0_alt_idx  = [222.35762166976929 266.0687934860357 -327.48699795456002 -383.16551638122462 89.743675980469447 -61.705085396679578 198.74117091566461 -329.9905668009938 -158.90935074534525 3.4979485816091369;-107.40944528579712 284.98732714159564 109.50484721664003 -353.20584214442096 -121.88649281109099 -146.82445100654729 4.8858705745241764 -115.32715990252382 344.97870245505123 123.68340401367681;308.72518682479858 325.52704211779553 -193.51869505535998 -323.24616790761729 -333.51666160265131 -284.32496468402593 -229.51432918437106 99.336246995946169 171.88738647018221 243.86885944574465;-21.041880130767822 272.37497137122239 238.36959571967998 -389.15745122858533 185.66238557664428 -369.44433029389364 265.89366057632105 -346.50313656241457 -1.2039295146868199 364.05431487781232;139.91585493087769 291.29350502678233 -64.653946552320008 -359.19777699178167 -31.504216669240805 226.39122897518041 72.038360235180733 -131.83972966394458 -174.29524549955585 -202.5342607305069;-189.85121202468872 310.21203868234227 -393.19526080511997 -342.93395383465969 -182.23361746322985 88.890715297701718 -121.81694010595976 82.823677234525405 329.59280770084069 -82.348805298439203;351.90896940231323 254.65761445569802 69.314356346880061 -312.97427959785603 -393.86378625479017 3.7713496878340607 -356.21713986485497 -397.77901108472116 156.50149171597167 37.836650133628495;22.141902446746826 295.19732943189786 -259.22695790592002 -378.88556291882406 125.31526092450542 -133.72916398964458 139.19084989583712 -183.1156041862512 -16.58982426889736 158.02210556569628;181.13673830032349 314.1158630874578 203.28265924608002 -348.9258886820204 -86.314907867055013 -218.84852959951232 -54.66445044530326 31.547802712218811 -219.11502587051709 278.20756099776395];   % 9 alternative starting points (Halton indices 1022 + 63 k)
PB.core_pidx   = 22;
PB.func_idx    = @(z) objective(z, D);
PB.func_F      = @(xr) objective(cell2idx(xr, PB), D);
PB.idx2cell    = @(z) idx2cell(z, PB);
PB.cell2idx    = @(xr) cell2idx(xr, PB);
PB.decode      = @(xr) decode_to_original(cell2idx(xr, PB), D);
fprintf('griewank_c_10D: n = 10, pattern C, m = 0, f* = 0 (E)\n');
end

function f = objective(z, D)
f = D.fn(decode_to_original(z, D));
end

function xh = decode_to_original(z, D)
if ~isnumeric(z) || numel(z) ~= 10 || ~all(isfinite(z)), error('griewank_c_10D: z must be a finite numeric vector with 10 entries.'); end
xh = zeros(10, 1);
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('griewank_c_10D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
v = z(2); tl = 1e-12 * (D.ubw(2) - D.lbw(2)); if v < D.lbw(2) - tl || v > D.ubw(2) + tl, error('griewank_c_10D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
v = min(max(v, D.lbw(2)), D.ubw(2));
if abs(v - D.cw(2)) <= 1e-12 * (D.ubw(2) - D.lbw(2)), xh(2) = D.cx(2); else, xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2)); end
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('griewank_c_10D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('griewank_c_10D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('griewank_c_10D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
v = z(6); tl = 1e-12 * (D.ubw(6) - D.lbw(6)); if v < D.lbw(6) - tl || v > D.ubw(6) + tl, error('griewank_c_10D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
v = min(max(v, D.lbw(6)), D.ubw(6));
if abs(v - D.cw(6)) <= 1e-12 * (D.ubw(6) - D.lbw(6)), xh(6) = D.cx(6); else, xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6)); end
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('griewank_c_10D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
v = z(8); tl = 1e-12 * (D.ubw(8) - D.lbw(8)); if v < D.lbw(8) - tl || v > D.ubw(8) + tl, error('griewank_c_10D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
v = min(max(v, D.lbw(8)), D.ubw(8));
if abs(v - D.cw(8)) <= 1e-12 * (D.ubw(8) - D.lbw(8)), xh(8) = D.cx(8); else, xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8)); end
v = z(9); tl = 1e-12 * (D.ubw(9) - D.lbw(9)); if v < D.lbw(9) - tl || v > D.ubw(9) + tl, error('griewank_c_10D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
v = min(max(v, D.lbw(9)), D.ubw(9));
if abs(v - D.cw(9)) <= 1e-12 * (D.ubw(9) - D.lbw(9)), xh(9) = D.cx(9); else, xh(9) = map_work_to_orig(v, D.lbw(9), D.ubw(9), D.lb(9), D.ub(9)); end
v = z(10); tl = 1e-12 * (D.ubw(10) - D.lbw(10)); if v < D.lbw(10) - tl || v > D.ubw(10) + tl, error('griewank_c_10D: x10 = %g outside [%g, %g].', v, D.lbw(10), D.ubw(10)); end
v = min(max(v, D.lbw(10)), D.ubw(10));
if abs(v - D.cw(10)) <= 1e-12 * (D.ubw(10) - D.lbw(10)), xh(10) = D.cx(10); else, xh(10) = map_work_to_orig(v, D.lbw(10), D.ubw(10), D.lb(10), D.ub(10)); end
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
xr = cell(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C', xr{i} = z(i);
        case 'D', xr{i} = PB.ProblemData{i}{1}(check_idx(z(i), PB.m, PB.name, i));
        case 'K', xr{i} = PB.ProblemData{i}{check_idx(z(i), PB.m, PB.name, i)};
    end
end
end

function z = cell2idx(xr, PB)
if ~iscell(xr) || numel(xr) ~= 10, error('griewank_c_10D: xr must be a cell array with 10 entries.'); end
z = zeros(1, 10);
for i = 1:10
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_c_10D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('griewank_c_10D: x%d must be a grid value.', i); end
            [d, k] = min(abs(grid - v)); if d > 1e-10 * max(1, max(abs(grid))), error('griewank_c_10D: x%d = %g is not a grid value.', i, v); end
            z(i) = k;
        case 'K'
            lab = xr{i}; if isstring(lab), lab = char(lab); end
            k = find(strcmp(lab, PB.ProblemData{i}), 1); if isempty(k), error('griewank_c_10D: x%d unknown label %s.', i, lab); end
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

