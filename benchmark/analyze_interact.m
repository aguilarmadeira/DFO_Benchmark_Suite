function T = analyze_interact(solver, folder, varargin)
%ANALYZE_INTERACT  Path statistics on folders with categorical variables (e.g. core/BK_MK16, core/BK_binary).
%
%   T = analyze_interact('my_solver', 'core/BK_MK16')
%   T = analyze_interact(solver, folder, 'results_dir', 'results', 'abs_levels', [], 'out_file', '')
%
% Requires runs made with run_benchmark(..., 'save_points', true). For every run:
%   evals_leave_x0cat   first evaluation whose categorical part differs from x0's
%                       (all categorical coordinates considered)
%   evals_first_optcat  first evaluation in the optimal category/combination (reference_optima.csv)
%   evals_to_<level>    first evaluation with f <= level, for ABSOLUTE levels common to the
%                       lambda-variants: b (floor of a wrong category; PB.b), b/10, 1e-6, 1e-9
%   evals_to_tau<τ>     Moré-Wild relative levels f <= tau f(x0), tau = 1e-1, 1e-3, 1e-5
%   effort per category: evals_in_optcat (evaluations whose categorical part is the optimal
%                       one), evals_wrong_after_first_opt (evaluations in other categories
%                       after the optimal one was first visited), evals_in_optcat_to_1e-6
%                       (evaluations spent inside the optimal category until f <= 1e-6: the
%                       continuous difficulty within the right category), n_categories_visited
% These are properties of the path of ANY solver on these problems (where it goes and
% when), not diagnostics of a particular mechanism. Output: table printed and written to
% <results_dir>/<solver>/<folder>/analyze_interact.csv, one row per run, with the
% instance parameters (lambda, obstacle ratio or regime) for pairing across variants.
opts = struct('results_dir', 'results', 'abs_levels', [], 'out_file', '');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
[~, fname] = fileparts(folder);                                % run_benchmark stores results under the folder's base name
d = fullfile(opts.results_dir, solver, fname);
L = dir(fullfile(d, '*_points.csv'));
if isempty(L), error('analyze_interact: no *_points.csv in %s (run with ''save_points'', true).', d); end
addpath(folder);
if isempty(opts.out_file), opts.out_file = fullfile(d, 'analyze_interact.csv'); end
fid = fopen(opts.out_file, 'w');
fprintf(fid, 'instance,seed,lambda,param,param_value,h_x0,f_x0,nevals,evals_leave_x0cat,evals_first_optcat,evals_to_b,evals_to_b_over_10,evals_to_1e-6,evals_to_1e-9,evals_to_tau1e-1,evals_to_tau1e-3,evals_to_tau1e-5,evals_in_optcat,evals_wrong_after_first_opt,evals_in_optcat_to_1e-6,n_categories_visited,fbest\n');
T = struct([]);
for k = 1:numel(L)
    tok = regexp(L(k).name, '^(.*)__s(\d+)_points\.csv$', 'tokens'); nm = tok{1}{1}; seed = str2double(tok{1}{2});
    try, [~, PB] = evalc('feval(nm)'); catch, evalc('PB = feval(nm);'); end
    PB = attach_reference_optimum(PB, folder);
    if ~isfield(PB, 'x_star_idx'), error('analyze_interact: no optimal solution for %s (PB.x_star_idx or reference_optima.csv).', nm); end
    P = dlmread(fullfile(d, L(k).name), ',', 1, 0); Z = P(:, 1:end-1); f = P(:, end); ne = size(P, 1);
    isK = strcmp(PB.var_types, 'K'); n = PB.nvar;
    if isfield(PB, 'x0_idx'), z0 = PB.x0_idx(:)';
    else
        z0 = zeros(1, n);
        for i = 1:n
            if isinf(PB.cardinality(i)), z0(i) = 0.5 * (PB.lb_C_work(i) + PB.ub_C_work(i)); else, z0(i) = round(0.5*(PB.cardinality(i)-1)) + 1; end
        end
    end
    f0 = f(1);
    leave = find(any(Z(:, isK) ~= repmat(z0(isK), ne, 1), 2), 1); if isempty(leave), leave = NaN; end
    optk = find(all(Z(:, isK) == repmat(PB.x_star_idx(isK), ne, 1), 2), 1); if isempty(optk), optk = NaN; end
    fb = cummin(f);
    b = NaN; if isfield(PB, 'b'), b = PB.b; end
    lev = [b, b/10, 1e-6, 1e-9]; e_abs = NaN(1, 4);
    for q = 1:4, j = find(fb <= lev(q), 1); if ~isempty(j), e_abs(q) = j; end, end
    e_rel = NaN(1, 3); taus = [1e-1 1e-3 1e-5];
    for q = 1:3, j = find(fb <= taus(q) * f0, 1); if ~isempty(j), e_rel(q) = j; end, end
    inopt = all(Z(:, isK) == repmat(PB.x_star_idx(isK), ne, 1), 2);
    n_in_opt = sum(inopt);
    if isnan(optk), n_wrong_after = NaN; else, n_wrong_after = sum(~inopt(optk:end)); end
    j6 = find(fb <= 1e-6, 1); if isempty(j6), n_opt_to6 = NaN; else, n_opt_to6 = sum(inopt(1:j6)); end
    [~, ~, cat_id] = unique(Z(:, isK), 'rows'); n_cats = max(cat_id);
    hx0 = sum(z0(isK) == PB.x_star_idx(isK));                 % categorical coordinates of x0 equal to the optimum's
    if isfield(PB, 'obstacle_ratio'), pname = 'obstacle_ratio'; pval = PB.obstacle_ratio;
    elseif isfield(PB, 'regime'), pname = 'regime'; pval = PB.regime;
    else, pname = 'none'; pval = NaN; end
    lam = NaN; if isfield(PB, 'lambda'), lam = PB.lambda; end
    fprintf(fid, '%s,%d,%.6g,%s,%.6g,%g,%.10g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%g,%g,%d,%.10g\n', nm, seed, lam, pname, pval, hx0, f0, ne, leave, optk, e_abs, e_rel, n_in_opt, n_wrong_after, n_opt_to6, n_cats, fb(end));
    T(end+1).instance = nm; T(end).seed = seed; T(end).lambda = lam; T(end).param = pval; T(end).f0 = f0; %#ok<AGROW>
    T(end).leave_x0cat = leave; T(end).first_optcat = optk; T(end).to_abs = e_abs; T(end).to_rel = e_rel; T(end).fbest = fb(end);
    T(end).in_optcat = n_in_opt; T(end).wrong_after_first = n_wrong_after; T(end).in_optcat_to_1e6 = n_opt_to6; T(end).n_cats = n_cats; T(end).h_x0 = hx0;
end
fclose(fid); rmpath(folder);
fprintf('%-36s %6s %7s %8s %9s %7s %8s %8s %8s %8s %6s\n', 'instance', 'lambda', 'param', 'leave_x0', 'first_opt', 'to_b', 'to_1e-6', 'in_opt', 'wrong>1st', 'opt_to6', 'ncats');
for r = 1:numel(T)
    fprintf('%-36s %6.3g %7.3g %8g %9g %7g %8g %8d %8g %8g %6d\n', T(r).instance, T(r).lambda, T(r).param, T(r).leave_x0cat, T(r).first_optcat, T(r).to_abs(1), T(r).to_abs(3), T(r).in_optcat, T(r).wrong_after_first, T(r).in_optcat_to_1e6, T(r).n_cats);
end
fprintf('written to %s\n', opts.out_file);
end
