function S = run_benchmark(solver, folder, varargin)
%RUN_BENCHMARK  Run a solver on every instance of a suite folder under the protocol.
%
%   S = run_benchmark(solver, folder)
%   S = run_benchmark(solver, folder, 'budget_per_n', 1000, 'seeds', 1:10, ...
%                     'results_dir', 'results', 'instances', {...}, 'verbose', 1)
%
% solver : struct with
%            .name  short identifier used in file names (e.g. 'my_solver')
%            .fun   handle  fun(P)  where P is the counted problem returned by
%                   bench_problem (P.func_F / P.func_idx, P.x0_idx, P.budget,
%                   P.seed, ...). The solver evaluates through P only; its
%                   return value is ignored (the trace is in the recorder).
%            .deterministic  (optional, default false) run once per instance,
%                   seed 1, whatever 'seeds' says.
% folder : a folder of the suite (e.g. core/D, core/CDK_baseline).
%
% Protocol (see PROTOCOL.md): budget = budget_per_n * n evaluations per run;
% x0 (centre) is evaluated first by this function, so every trace starts with
% (1, f(x0)); seeds 1:10 for stochastic solvers; one
% file per run  results/<solver>/<folder>/<instance>__s<seed>.csv  with the
% trace of improvements (neval, fbest), and one summary per (solver, folder)
% results/<solver>/<folder>/summary.csv with one row per run:
%   instance,n,seed,budget,nevals,n_distinct,fbest,f_star,gap,n_nonfinite,time_s
% (nevals counts calls to the objective, n_distinct the distinct points: a
% solver without a cache may re-evaluate x0 or other points)
% Existing run files are skipped (delete them to rerun). With 'save_points', true
% the whole sequence of evaluated points is written to <instance>__s<seed>_points.csv
% (index representation), for analyses of the path (e.g. analyze_interact).
opts = struct('budget_per_n', 1000, 'seeds', 1:10, 'results_dir', 'results', 'instances', {{}}, 'verbose', 1, 'save_points', false);
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
if ~isfield(solver, 'deterministic'), solver.deterministic = false; end
if solver.deterministic, seeds = 1; else, seeds = opts.seeds; end

addpath(folder);
[~, fname] = fileparts(folder);
out = fullfile(opts.results_dir, solver.name, fname);
if ~exist(out, 'dir'), mkdir(out); end
L = dir(fullfile(folder, '*.m'));
L = L(~ismember({L.name}, {'add_problems.m', 'add_strategy.m', 'add_all_strategies.m'}));   % every wrapper of the folder
names = cellfun(@(s) s(1:end-2), {L.name}, 'UniformOutput', false);
if ~isempty(opts.instances), names = intersect(names, opts.instances, 'stable'); end

sumfile = fullfile(out, 'summary.csv');
if ~exist(sumfile, 'file')
    fid = fopen(sumfile, 'w'); fprintf(fid, 'instance,n,seed,budget,nevals,n_distinct,fbest,f_star,gap,n_nonfinite,time_s\n'); fclose(fid);
end
S = struct([]);
for k = 1:numel(names)
    PB = load_quiet(names{k});
    n = PB.nvar; budget = opts.budget_per_n * n;
    for seed = seeds
        runfile = fullfile(out, sprintf('%s__s%d.csv', names{k}, seed));
        if exist(runfile, 'file'), continue; end
        [P, rec] = bench_problem(PB, budget, seed, opts.save_points);
        t0 = tic;
        f0 = P.func_F(P.x0);                 % protocol: x0 is always the first evaluation
        if ~isfinite(f0), fprintf('  WARNING %s: f(x0) is not finite; the convergence test of this instance is undefined\n', names{k}); end
        try
            solver.fun(P);
        catch ME
            if ~strcmp(ME.identifier, 'BENCH:budget')
                fprintf('  ERROR in %s on %s (seed %d): %s\n', solver.name, names{k}, seed, ME.message);
                continue;
            end
        end
        t = toc(t0);
        if isfinite(PB.f_star), gap = rec.fbest - PB.f_star; else, gap = NaN; end
        fid = fopen(runfile, 'w');
        fprintf(fid, '# instance=%s n=%d seed=%d budget=%d f_star=%.16g f0=%.16g solver=%s\n', names{k}, n, seed, budget, PB.f_star, f0, solver.name);
        fprintf(fid, 'neval,fbest\n');
        fprintf(fid, '%d,%.16g\n', rec.trace');
        fclose(fid);
        if opts.save_points && ~isempty(rec.points)
            pts = rec.points(1:rec.nevals, :);
            pfile = fullfile(out, sprintf('%s__s%d_points.csv', names{k}, seed));
            fid = fopen(pfile, 'w'); fprintf(fid, '# instance=%s seed=%d columns: z(1..n) then f (index representation; continuous coordinates as values)\n', names{k}, seed);
            fmt = [repmat('%.17g,', 1, size(pts, 2) - 1) '%.17g\n']; fprintf(fid, fmt, pts'); fclose(fid);
        end
        fid = fopen(sumfile, 'a');
        fprintf(fid, '%s,%d,%d,%d,%d,%d,%.16g,%.16g,%.6g,%d,%.2f\n', names{k}, n, seed, budget, rec.nevals, rec.n_distinct, rec.fbest, PB.f_star, gap, rec.n_nonfinite, t);
        fclose(fid);
        if opts.verbose
            fprintf('%-45s n=%2d seed %2d: f=%.10g f*=%.10g gap=%.3g  %d evals (%d distinct)  %.1fs\n', names{k}, n, seed, rec.fbest, PB.f_star, gap, rec.nevals, rec.n_distinct, t);
        end
        S(end+1).instance = names{k}; S(end).n = n; S(end).seed = seed; S(end).fbest = rec.fbest; %#ok<AGROW>
        S(end).f_star = PB.f_star; S(end).gap = gap; S(end).nevals = rec.nevals; S(end).n_distinct = rec.n_distinct; S(end).trace = rec.trace; S(end).time = t;
    end
end
rmpath(folder);
end

function PB = load_quiet(name)
try
    [~, PB] = evalc('feval(name)');     % MATLAB
catch
    evalc('PB = feval(name);');         % Octave
end
end
