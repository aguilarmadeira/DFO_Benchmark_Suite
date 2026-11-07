function R = data_profiles(solvers, folder, varargin)
%DATA_PROFILES  Data and performance profiles (More & Wild, 2009) from run_benchmark results.
%
%   R = data_profiles({'random_search','my_solver'}, 'core/D')
%   R = data_profiles(solvers, folder, 'taus', [1e-1 1e-3 1e-5 1e-7], ...
%                     'results_dir', 'results', 'out_dir', 'profiles', 'plot', true)
%
% Problems: every (instance, seed) pair present for ALL listed solvers (a
% deterministic solver's single run stands for every seed). Convergence test
% of More & Wild: f <= f_L + tau * (f0 - f_L), where f0 = f(x0) (first
% evaluation of every run) and f_L the TARGET of the instance: f_star when
% known, otherwise a best-known value. Targets are read from a versioned file
% targets/<folder>_targets.csv (option 'targets'), so that adding a solver to
% a comparison never changes the targets of earlier comparisons. When no
% targets file exists the best value found by the listed solvers is used and
% the table is written to <out_dir>/<folder>_targets_auto.csv: review it, copy
% it to targets/ and version it (see PROTOCOL.md). Instances flagged trivial
% in <folder>/INDEX.csv (x0 is an optimum) are excluded unless
% 'exclude_trivial' is false.
% Data profile: d_s(alpha) = fraction of problems solved within alpha*(n+1)
% evaluations. Performance profile at the full budget: rho_s(r) = fraction
% of problems where t_s <= r * min_s t_s (unsolved: t = Inf).
% Outputs (out_dir): <folder>_tau<tau>_data.csv, <folder>_tau<tau>_perf.csv,
% <folder>_instances.csv (per instance and solver: fraction of seeds solved
% at each tau and median evaluations), and PNG plots when 'plot' is true.
opts = struct('taus', [1e-1 1e-3 1e-5 1e-7], 'results_dir', 'results', 'out_dir', 'profiles', 'plot', true, 'alpha_max', [], 'targets', '', 'exclude_trivial', true, 'index_file', '', 'instances', {{}});
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
if ischar(solvers), solvers = {solvers}; end
[~, fname] = fileparts(folder);
if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end

% ---- read all runs ----
runs = struct('solver', {}, 'instance', {}, 'seed', {}, 'n', {}, 'budget', {}, 'f_star', {}, 'f0', {}, 'trace', {});
for s = 1:numel(solvers)
    d = fullfile(opts.results_dir, solvers{s}, fname);
    L = dir(fullfile(d, '*__s*.csv')); L = L(cellfun(@isempty, regexp({L.name}, '_points\.csv$', 'once')));   % traces only (not the save_points files)
    if isempty(L), error('data_profiles: no runs in %s', d); end
    for k = 1:numel(L)
        [hdr, tr] = read_run(fullfile(d, L(k).name));
        runs(end+1) = struct('solver', solvers{s}, 'instance', hdr.instance, 'seed', hdr.seed, 'n', hdr.n, ...
            'budget', hdr.budget, 'f_star', hdr.f_star, 'f0', hdr.f0, 'trace', tr); %#ok<AGROW>
    end
end
instances = unique({runs.instance});
if ~isempty(opts.instances)                       % restrict to a given list (e.g. paired_nontrivial)
    instances = intersect(instances, opts.instances, 'stable'); runs = runs(ismember({runs.instance}, instances));
end
% trivial instances (x0 is an optimum, declared in <folder>/INDEX.csv by validate_suite) are excluded by default
if opts.exclude_trivial
    if isempty(opts.index_file), opts.index_file = fullfile(folder, 'INDEX.csv'); end
    triv = read_trivial(opts.index_file);
    if ~isempty(triv)
        drop = ismember(instances, triv);
        if any(drop), fprintf('data_profiles: %d trivial instance(s) excluded (x0 is an optimum): %s\n', sum(drop), strjoin(instances(drop), ', ')); end
        instances = instances(~drop); runs = runs(ismember({runs.instance}, instances));
    end
end
seeds = unique([runs.seed]);

% ---- f_L per instance: versioned targets, else pooled best (written for review) ----
if isempty(opts.targets), opts.targets = fullfile('targets', [fname '_targets.csv']); end
tgt = read_targets(opts.targets);
fL = struct(); src = struct(); n_auto = 0;
for i = 1:numel(instances)
    rr = runs(strcmp({runs.instance}, instances{i})); k = key(instances{i});
    if isfield(tgt, k)
        fL.(k) = tgt.(k); src.(k) = 'targets_file';
    elseif isfinite(rr(1).f_star)
        fL.(k) = rr(1).f_star; src.(k) = 'f_star';
    else
        fL.(k) = min(cellfun(@(t) min([t(:,2); Inf]), {rr.trace})); src.(k) = 'best_of_listed_solvers'; n_auto = n_auto + 1;
    end
end
if n_auto > 0
    afile = fullfile(opts.out_dir, [fname '_targets_auto.csv']);
    fid = fopen(afile, 'w'); fprintf(fid, 'instance,f_target,source\n');
    for i = 1:numel(instances), k = key(instances{i}); fprintf(fid, '%s,%.16g,%s\n', instances{i}, fL.(k), src.(k)); end
    fclose(fid);
    fprintf('data_profiles: %d instance(s) without f_star and without a versioned target use the best value of the listed solvers.\n  Review %s, copy it to %s and version it so that later comparisons use the same targets.\n', n_auto, afile, opts.targets);
end

% ---- t(problem, solver, tau): evaluations to convergence ----
nP = numel(instances) * numel(seeds); nS = numel(solvers); nT = numel(opts.taus);
T = Inf(nP, nS, nT); Nn = zeros(nP, 1); B = zeros(nP, 1); p = 0; labels = cell(nP, 1);
for i = 1:numel(instances)
    for sd = seeds
        p = p + 1; labels{p} = sprintf('%s__s%d', instances{i}, sd);
        for s = 1:nS
            r = runs(strcmp({runs.solver}, solvers{s}) & strcmp({runs.instance}, instances{i}) & [runs.seed] == sd);
            if isempty(r)   % deterministic solver: its single run stands for every seed
                r = runs(strcmp({runs.solver}, solvers{s}) & strcmp({runs.instance}, instances{i}));
                if isempty(r), error('data_profiles: %s has no run on %s', solvers{s}, instances{i}); end
                r = r(1);
            end
            Nn(p) = r.n; B(p) = r.budget; fl = fL.(key(instances{i}));
            for t = 1:nT
                thr = fl + opts.taus(t) * (r.f0 - fl);
                j = find(r.trace(:,2) <= thr, 1);
                if ~isempty(j), T(p, s, t) = r.trace(j, 1); end
            end
        end
    end
end

% ---- data profiles ----
if isempty(opts.alpha_max), opts.alpha_max = max(B ./ (Nn + 1)); end
alpha = unique([0, logspace(0, log10(opts.alpha_max), 200)]);
R = struct('solvers', {solvers}, 'taus', opts.taus, 'alpha', alpha, 'problems', {labels}, 'T', T);
R.data = zeros(numel(alpha), nS, nT); R.solved = zeros(nS, nT);
for t = 1:nT
    for s = 1:nS
        tt = T(:, s, t) ./ (Nn + 1);
        R.data(:, s, t) = arrayfun(@(a) mean(tt <= a), alpha);
        R.solved(s, t) = mean(isfinite(T(:, s, t)));
    end
    write_table(fullfile(opts.out_dir, sprintf('%s_tau%g_data.csv', fname, opts.taus(t))), 'alpha', alpha(:), solvers, squeeze(R.data(:, :, t)));
end
% ---- performance profiles at full budget ----
rgrid = unique([1, logspace(0, log10(max(B)), 200)]);
R.rgrid = rgrid; R.perf = zeros(numel(rgrid), nS, nT);
for t = 1:nT
    tmin = min(T(:, :, t), [], 2);
    for s = 1:nS
        ratio = T(:, s, t) ./ tmin; ratio(~isfinite(tmin)) = Inf;
        R.perf(:, s, t) = arrayfun(@(r) mean(ratio <= r), rgrid);
    end
    write_table(fullfile(opts.out_dir, sprintf('%s_tau%g_perf.csv', fname, opts.taus(t))), 'ratio', rgrid(:), solvers, squeeze(R.perf(:, :, t)));
end
% ---- per-instance table ----
fid = fopen(fullfile(opts.out_dir, sprintf('%s_instances.csv', fname)), 'w');
fprintf(fid, 'instance,n,solver');
for t = 1:nT, fprintf(fid, ',solved_tau%g,median_evals_tau%g', opts.taus(t), opts.taus(t)); end
fprintf(fid, '\n');
for i = 1:numel(instances)
    rows = (i-1)*numel(seeds) + (1:numel(seeds));
    for s = 1:nS
        fprintf(fid, '%s,%d,%s', instances{i}, Nn(rows(1)), solvers{s});
        for t = 1:nT
            tt = T(rows, s, t); fprintf(fid, ',%.3f,%g', mean(isfinite(tt)), median(tt));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
% ---- console summary ----
fprintf('\n%s: %d instances x %d seeds = %d problems; fraction solved at the full budget\n', fname, numel(instances), numel(seeds), nP);
fprintf('%-28s', 'solver'); fprintf('  tau=%-7g', opts.taus); fprintf('\n');
for s = 1:nS, fprintf('%-28s', solvers{s}); fprintf('  %-11.3f', R.solved(s, :)); fprintf('\n'); end
% ---- plots ----
if opts.plot
    try
        for t = 1:nT
            fh = figure('visible', 'off'); hold on;
            for s = 1:nS, semilogx(alpha(2:end), R.data(2:end, s, t), 'LineWidth', 1.5); end
            xlabel('\alpha  (evaluations / (n+1))'); ylabel('fraction of problems solved'); ylim([0 1]); grid on;
            title(sprintf('%s, data profile, \\tau = %g', strrep(fname, '_', '\_'), opts.taus(t)));
            legend(strrep(solvers, '_', '\_'), 'Location', 'southeast');
            print(fh, fullfile(opts.out_dir, sprintf('%s_tau%g_data.png', fname, opts.taus(t))), '-dpng', '-r120'); close(fh);
        end
    catch ME
        fprintf('(plots skipped: %s)\n', ME.message);
    end
end
end

function [hdr, tr] = read_run(file)
fid = fopen(file, 'r'); line = fgetl(fid);
hdr = struct('instance', '', 'seed', NaN, 'n', NaN, 'budget', NaN, 'f_star', NaN, 'f0', NaN);
toks = regexp(line, '(\w+)=(\S+)', 'tokens');
for k = 1:numel(toks)
    switch toks{k}{1}
        case 'instance', hdr.instance = toks{k}{2};
        case {'seed','n','budget','f_star','f0'}, hdr.(toks{k}{1}) = str2double(toks{k}{2});
    end
end
fgetl(fid);                               % column header
tr = fscanf(fid, '%f,%f', [2 Inf])'; fclose(fid);
if isempty(tr), tr = zeros(0, 2); end
end

function triv = read_trivial(file)
triv = {};
if exist(file, 'file') ~= 2, return; end
fid = fopen(file, 'r'); hdr = strsplit(strtrim(fgetl(fid)), ','); col = find(strcmp(hdr, 'x0_is_optimum'), 1);
if isempty(col), fclose(fid); return; end
while true
    line = fgetl(fid); if ~ischar(line), break; end
    c = strsplit(line, ','); if numel(c) >= col && strcmp(strtrim(c{col}), '1'), triv{end+1} = c{1}; end %#ok<AGROW>
end
fclose(fid);
end

function tgt = read_targets(file)
tgt = struct();
if exist(file, 'file') ~= 2, return; end
fid = fopen(file, 'r'); fgetl(fid);
while true
    line = fgetl(fid); if ~ischar(line), break; end
    c = strsplit(line, ','); if numel(c) < 2, continue; end
    tgt.(key(c{1})) = str2double(c{2});
end
fclose(fid);
end

function k = key(name)
k = ['i_' regexprep(name, '\W', '_')];
if numel(k) > 63, k = k(1:63); end
end

function write_table(file, xname, x, solvers, Y)
fid = fopen(file, 'w'); fprintf(fid, '%s', xname); fprintf(fid, ',%s', solvers{:}); fprintf(fid, '\n');
for i = 1:numel(x), fprintf(fid, '%.6g', x(i)); fprintf(fid, ',%.4f', Y(i, :)); fprintf(fid, '\n'); end
fclose(fid);
end
