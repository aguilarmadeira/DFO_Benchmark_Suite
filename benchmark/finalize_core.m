function finalize_core(core_root, varargin)
%FINALIZE_CORE  Validate every folder of the version-2 core and write its versioned targets.
%
%   finalize_core                       % ./core, enumerate_max 2e6
%   finalize_core('core', 'enumerate_max', 2e6, 'n_random', 200, 'targets_dir', '<core_root>/../targets')
%
% For each sub-folder of core_root: validate_suite (writes INDEX.csv; checks admissibility of x0,
% the alternative starting points and a Halton sample, F(x*) = f* through reference_optima.csv,
% local minimality of x* in index space, enumeration of the finite domain when it has at most
% enumerate_max points, trivial instances) and make_targets (targets/<folder>_targets.csv).
% Prints one summary line per folder and a final table. Run from the kit root with
% addpath('benchmark'); the wrappers are autonomous.
if nargin < 1 || isempty(core_root), core_root = 'core'; end
opts = struct('enumerate_max', 2e6, 'n_random', 200, 'targets_dir', '');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
d = dir(core_root); if isempty(d), error('finalize_core: folder "%s" not found (current folder %s).', core_root, pwd); end
core_root = d(1).folder;                                                     % absolute path
if isempty(opts.targets_dir), opts.targets_dir = fullfile(fileparts(core_root), 'targets'); end   % next to core/
L = dir(core_root); L = L([L.isdir]); L = L(~ismember({L.name}, {'.', '..'}));
if isempty(L)
    error('finalize_core: no sub-folders found in "%s" (current folder: %s). Pass the path of the core/ folder, e.g. finalize_core(''C:\\...\\DFO_Benchmark_Suite_v2.0.0\\core'').', core_root, pwd);
end
S = cell(numel(L), 6);
t_all = tic;
for k = 1:numel(L)
    folder = fullfile(core_root, L(k).name);
    t = tic;
    V = validate_suite(folder, 'enumerate_max', opts.enumerate_max, 'n_random', opts.n_random);
    make_targets(folder, 'targets_dir', opts.targets_dir);
    npb = sum(arrayfun(@(v) ~isempty(v.problems), V));
    S(k, :) = {L(k).name, numel(V), npb, sum([V.certified]), sum([V.x0_is_optimum]), toc(t)};
end
fprintf('\n%-16s %10s %10s %10s %10s %10s\n', 'folder', 'instances', 'problems', 'certified', 'trivial', 'seconds');
for k = 1:size(S, 1), fprintf('%-16s %10d %10d %10d %10d %10.0f\n', S{k, :}); end
fprintf('total %.0f s; INDEX.csv in each folder, targets in %s. Publish only if every "problems" and "trivial" is 0.\n', toc(t_all), opts.targets_dir);
end
