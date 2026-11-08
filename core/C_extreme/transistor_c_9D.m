function PB = transistor_c_9D()
%TRANSISTOR_C_9D  DFO Benchmark Suite v2.1.0 core: transistor, dimension 9, pattern C, strategy extreme.
%
% Copyright (C) 2026 J. F. A. Madeira
%
% Objective F(z) = f(T(z)), where f is the original function (embedded below) and T acts on the finite
% coordinates only: uniform grid of 0 levels on the interval of the domain, piecewise-affine map with |T - Id| <=
% half a cell, identity on the continuous coordinates (up to the affine scaling of the strategy).
% Global optimum of F on this domain: f* = 1.60049326654736e-24 (E: exact by structure); the optimal solution is listed in
% reference_optima.csv of this folder, not here. Starting point of the protocol: PB.x0_idx (Halton).
% Variables: x1:C x2:C x3:C x4:C x5:C x6:C x7:C x8:C x9:C
%
lb_orig = [-10 -10 -10 -10 -7.3500853766819212 -6.9926613213823039 -10 -10 -10];
ub_orig = [6.8211669921875 7.0331250317532898 6.7027327999999997 6.1800083298625577 10 10 6.7876079069934505 6.4033133570184386 7.3164625626695159];
lb_C_work = [-10 -10 -10 -10 -735008537.66819215 -699266132.13823032 -1000000000 -1000000000 -999999999.99999988];
ub_C_work = [6.8211669921875 7.0331250317532898 6.7027327999999997 6.1800083298625585 1000000000 1000000000 678760790.69934511 640331335.70184398 731646256.26695144];
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
PB.suite_release = '2.1.0';
PB.description = 'DFO Benchmark Suite v2.1.0 core: C adaptation of transistor (F = f o T).';
PB.nvar        = 9;
PB.nobj        = 1;
PB.var_types   = {'C','C','C','C','C','C','C','C','C'};
PB.type_counts = struct('C', 9, 'D', 0, 'K', 0);
PB.proportion  = [9 0 0];
PB.pattern     = 'C';
PB.strategy    = 'extreme';
PB.scale_factors_C = [1 1 1 1 100000000 100000000 100000000 100000000 100000000];
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
D.cw = [0.89999995261685584 0.44998747198158462 1.0000064824653876 2.0000685416236834 799997144.0507853 799969268.42185426 500003127.59296727 99998772.345679045 200005248.3486191];   % anchors of the continuous maps (NaN for finite)
D.cx = [0.8999999526168565 0.44998747198158512 1.0000064824653869 2.000068541623683 7.9999714405078528 7.9996926842185427 5.0000312759296719 0.99998772345678977 2.0000524834861899];
D.fn = @(xh) transistor_orig(xh);
PB.f_star      = 1.6004932665473596e-24;
PB.f_star_justification = 'E';   % exact by structure
PB.n_minimizers_documented = 1;
PB.x0_idx      = [-3.4210181832313538 -9.2678949460517561 -1.6299265392639999 2.8981824003985572 336499364.63327587 511954060.36448646 -145072359.38738823 545388699.70918703 -698559482.96429658];   % protocol starting point: Halton point 1062 (index representation)
PB.x0_alt_idx  = [0.91568893194198608 -8.2164674749558735 -9.8877576355840002 3.6057629396053752 -207075812.08170241 300802647.49523902 -615590089.65260267 -644861928.33981073 587899834.07334387;-7.5606022477149963 -9.4080852755312065 1.7373443932160004 2.0490857533503757 -594226477.51179492 -40288096.370468736 488435172.8651228 -113230996.82425463 145986328.22071946;3.084042489528656 -9.1043395616590637 -6.6487636910079999 2.7566662925571919 582867969.90697122 -251439509.23971629 17917442.599908352 418399934.69130135 -295927177.63190532;-5.490810215473175 -8.6136734084809845 4.976338337791999 3.4642468317640116 195717304.47687876 -462590922.10896385 -452600287.66530621 -685756615.37946892 -737840683.48452973;-1.1541031002998352 -9.8052912090563176 -3.2814927585280005 1.5841042561573229 -362078282.35490745 885529636.97930908 -923118017.93052053 -154125683.86391282 548618633.55311072;-9.6303942799568176 -8.7538637379604367 -7.6910142177280001 2.2916847953641408 971440676.34874439 674378224.11006141 279316181.6361382 377505247.65164328 106705127.70048606;6.5336958765983582 -8.5124248371902702 0.085778173952000714 2.999265334570957 427865499.63376594 463226811.24081397 -289610485.67800975 -812745380.39735472 -335208378.15213847;-2.0904375910758972 -9.7040426377656033 -4.3237432852480007 3.706845873777775 40714834.203673482 122136067.37510633 -760128215.94322407 -281114448.88179874 -777121884.00476313;2.2462695240974426 -9.2133764845875259 3.453049106432001 2.150168687522779 -502860342.51130474 -89015345.494141221 442305983.62343478 250516482.63375735 434190788.5593878];   % 9 alternative starting points (Halton indices 1062 + 63 k)
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
v = z(1); tl = 1e-12 * (D.ubw(1) - D.lbw(1)); if v < D.lbw(1) - tl || v > D.ubw(1) + tl, error('transistor_c_9D: x1 = %g outside [%g, %g].', v, D.lbw(1), D.ubw(1)); end
v = min(max(v, D.lbw(1)), D.ubw(1));
if abs(v - D.cw(1)) <= 1e-12 * (D.ubw(1) - D.lbw(1)), xh(1) = D.cx(1); else, xh(1) = map_work_to_orig(v, D.lbw(1), D.ubw(1), D.lb(1), D.ub(1)); end
v = z(2); tl = 1e-12 * (D.ubw(2) - D.lbw(2)); if v < D.lbw(2) - tl || v > D.ubw(2) + tl, error('transistor_c_9D: x2 = %g outside [%g, %g].', v, D.lbw(2), D.ubw(2)); end
v = min(max(v, D.lbw(2)), D.ubw(2));
if abs(v - D.cw(2)) <= 1e-12 * (D.ubw(2) - D.lbw(2)), xh(2) = D.cx(2); else, xh(2) = map_work_to_orig(v, D.lbw(2), D.ubw(2), D.lb(2), D.ub(2)); end
v = z(3); tl = 1e-12 * (D.ubw(3) - D.lbw(3)); if v < D.lbw(3) - tl || v > D.ubw(3) + tl, error('transistor_c_9D: x3 = %g outside [%g, %g].', v, D.lbw(3), D.ubw(3)); end
v = min(max(v, D.lbw(3)), D.ubw(3));
if abs(v - D.cw(3)) <= 1e-12 * (D.ubw(3) - D.lbw(3)), xh(3) = D.cx(3); else, xh(3) = map_work_to_orig(v, D.lbw(3), D.ubw(3), D.lb(3), D.ub(3)); end
v = z(4); tl = 1e-12 * (D.ubw(4) - D.lbw(4)); if v < D.lbw(4) - tl || v > D.ubw(4) + tl, error('transistor_c_9D: x4 = %g outside [%g, %g].', v, D.lbw(4), D.ubw(4)); end
v = min(max(v, D.lbw(4)), D.ubw(4));
if abs(v - D.cw(4)) <= 1e-12 * (D.ubw(4) - D.lbw(4)), xh(4) = D.cx(4); else, xh(4) = map_work_to_orig(v, D.lbw(4), D.ubw(4), D.lb(4), D.ub(4)); end
v = z(5); tl = 1e-12 * (D.ubw(5) - D.lbw(5)); if v < D.lbw(5) - tl || v > D.ubw(5) + tl, error('transistor_c_9D: x5 = %g outside [%g, %g].', v, D.lbw(5), D.ubw(5)); end
v = min(max(v, D.lbw(5)), D.ubw(5));
if abs(v - D.cw(5)) <= 1e-12 * (D.ubw(5) - D.lbw(5)), xh(5) = D.cx(5); else, xh(5) = map_work_to_orig(v, D.lbw(5), D.ubw(5), D.lb(5), D.ub(5)); end
v = z(6); tl = 1e-12 * (D.ubw(6) - D.lbw(6)); if v < D.lbw(6) - tl || v > D.ubw(6) + tl, error('transistor_c_9D: x6 = %g outside [%g, %g].', v, D.lbw(6), D.ubw(6)); end
v = min(max(v, D.lbw(6)), D.ubw(6));
if abs(v - D.cw(6)) <= 1e-12 * (D.ubw(6) - D.lbw(6)), xh(6) = D.cx(6); else, xh(6) = map_work_to_orig(v, D.lbw(6), D.ubw(6), D.lb(6), D.ub(6)); end
v = z(7); tl = 1e-12 * (D.ubw(7) - D.lbw(7)); if v < D.lbw(7) - tl || v > D.ubw(7) + tl, error('transistor_c_9D: x7 = %g outside [%g, %g].', v, D.lbw(7), D.ubw(7)); end
v = min(max(v, D.lbw(7)), D.ubw(7));
if abs(v - D.cw(7)) <= 1e-12 * (D.ubw(7) - D.lbw(7)), xh(7) = D.cx(7); else, xh(7) = map_work_to_orig(v, D.lbw(7), D.ubw(7), D.lb(7), D.ub(7)); end
v = z(8); tl = 1e-12 * (D.ubw(8) - D.lbw(8)); if v < D.lbw(8) - tl || v > D.ubw(8) + tl, error('transistor_c_9D: x8 = %g outside [%g, %g].', v, D.lbw(8), D.ubw(8)); end
v = min(max(v, D.lbw(8)), D.ubw(8));
if abs(v - D.cw(8)) <= 1e-12 * (D.ubw(8) - D.lbw(8)), xh(8) = D.cx(8); else, xh(8) = map_work_to_orig(v, D.lbw(8), D.ubw(8), D.lb(8), D.ub(8)); end
v = z(9); tl = 1e-12 * (D.ubw(9) - D.lbw(9)); if v < D.lbw(9) - tl || v > D.ubw(9) + tl, error('transistor_c_9D: x9 = %g outside [%g, %g].', v, D.lbw(9), D.ubw(9)); end
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
if ~iscell(xr) || numel(xr) ~= 9, error('transistor_c_9D: xr must be a cell array with 9 entries.'); end
z = zeros(1, 9);
for i = 1:9
    switch PB.var_types{i}
        case 'C'
            v = xr{i}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('transistor_c_9D: x%d must be a finite scalar.', i); end
            z(i) = double(v);
        case 'D'
            v = xr{i}; grid = PB.ProblemData{i}{1}; if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), error('transistor_c_9D: x%d must be a grid value.', i); end
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

