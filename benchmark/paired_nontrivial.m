function [bases, per_folder] = paired_nontrivial(folders)
%PAIRED_NONTRIVIAL  Base functions that are non-trivial in EVERY listed folder.
%
%   [bases, per_folder] = paired_nontrivial({'core/CD_baseline', 'core/CD_moderate', 'core/CD_extreme'})
%
% Reads <folder>/INDEX.csv (validate_suite) and returns the base functions
% (column base_function) whose instance is not trivial (x0_is_optimum = 0) in
% all folders, and per_folder{i} = the instance names of folder i for those
% bases. Use the latter as the 'instances' option of data_profiles so that a
% comparison across configurations (proportions, cardinalities, strategies)
% is PAIRED: the same base functions in every configuration, and the effect
% of the configuration is not confounded with a change of the function set.
if ischar(folders), folders = {folders}; end
sets = cell(size(folders)); names = cell(size(folders)); bf = cell(size(folders));
for i = 1:numel(folders)
    T = read_index(fullfile(folders{i}, 'INDEX.csv'));
    ok = ~strcmp(T.x0_is_optimum, '1');
    sets{i} = T.base_function(ok); names{i} = T.instance(ok); bf{i} = T.base_function;
end
bases = sets{1};
for i = 2:numel(folders), bases = intersect(bases, sets{i}, 'stable'); end
per_folder = cell(size(folders));
for i = 1:numel(folders)
    [tf, loc] = ismember(bases, bf{i}); T = read_index(fullfile(folders{i}, 'INDEX.csv'));
    per_folder{i} = T.instance(loc(tf));
end
fprintf('paired_nontrivial: %d base function(s) non-trivial in all %d folders\n', numel(bases), numel(folders));
end

function S = read_index(file)
fid = fopen(file, 'r'); hdr = strsplit(strtrim(fgetl(fid)), ','); S = struct();
for k = 1:numel(hdr), S.(hdr{k}) = {}; end
while true
    line = fgetl(fid); if ~ischar(line), break; end
    c = strsplit(line, ',');
    for k = 1:numel(hdr), if k <= numel(c), S.(hdr{k}){end+1} = c{k}; else, S.(hdr{k}){end+1} = ''; end, end
end
fclose(fid);
if ~isfield(S, 'base_function'), error('paired_nontrivial: %s has no base_function column (rerun validate_suite).', file); end
end
