function T = make_targets(folder, varargin)
%MAKE_TARGETS  Build (or update) the versioned targets file of a folder from its INDEX.csv.
%
%   make_targets('core/D')
%   make_targets(folder, 'targets_dir', 'targets', 'best_known', 'profiles/D_targets_auto.csv')
%
% For every instance of <folder>/INDEX.csv (written by validate_suite):
%   source 'certified'   f_star_certified, when the optimum of the finite domain was
%                        certified by enumeration (also for instances without
%                        documented minimizer: e.g. cosine_mixture 2D, shubert 2D);
%   source 'f_star'      f at the injected documented minimizer otherwise, when known;
%   source 'best_known'  a value taken from the optional best-known file (instance,
%                        f_target, ...; e.g. the *_targets_auto.csv written by
%                        data_profiles), when neither is available;
%   source 'none'        NaN: the instance has no target yet (data_profiles then
%                        falls back to the best value of the listed solvers).
% An existing targets file is read first and its 'best_known' values are kept
% unless the new file gives a lower one; certified / f_star values always win.
% Writes targets/<folder>_targets.csv with a date column; version it with the results.
opts = struct('targets_dir', 'targets', 'best_known', '');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
[~, fname] = fileparts(folder);
tdir = opts.targets_dir;
if ~(ispc && ~isempty(regexp(tdir, '^[A-Za-z]:[\\/]', 'once')) || strncmp(tdir, '/', 1)), tdir = fullfile(pwd, tdir); end   % absolute: exist(name,'dir') would search the whole path
if ~exist(tdir, 'dir'), [ok, msg] = mkdir(tdir); if ~ok, error('make_targets: cannot create %s (%s)', tdir, msg); end, end
out = fullfile(tdir, [fname '_targets.csv']);
idx = read_csv(fullfile(folder, 'INDEX.csv'));
old = struct(); if exist(out, 'file') == 2, old = read_csv(out); end
bk  = struct(); if ~isempty(opts.best_known) && exist(opts.best_known, 'file') == 2, bk = read_csv(opts.best_known); end
today = datestr(now, 'yyyy-mm-dd');
T = struct('instance', {}, 'f_target', {}, 'source', {}, 'date', {});
for i = 1:numel(idx.instance)
    nm = idx.instance{i}; f = NaN; src = 'none'; dt = today;
    if str2double(idx.certified{i}) == 1 && isfinite(str2double(idx.f_star_certified{i}))
        f = str2double(idx.f_star_certified{i}); src = 'certified';
    elseif isfinite(str2double(idx.f_star{i}))
        f = str2double(idx.f_star{i}); src = 'f_star';
    else
        cand = NaN;
        if isfield(bk, 'instance'), j = find(strcmp(bk.instance, nm), 1); if ~isempty(j), cand = str2double(bk.f_target{j}); end, end
        if isfield(old, 'instance')
            j = find(strcmp(old.instance, nm), 1);
            if ~isempty(j) && strcmp(old.source{j}, 'best_known')
                fo = str2double(old.f_target{j});
                if ~(isfinite(cand) && cand < fo - 1e-12*max(1,abs(fo))), cand = fo; dt = old.date{j}; end
            end
        end
        if isfinite(cand), f = cand; src = 'best_known'; end
    end
    T(end+1) = struct('instance', nm, 'f_target', f, 'source', src, 'date', dt); %#ok<AGROW>
end
fid = fopen(out, 'w'); if fid < 0, error('make_targets: cannot write %s', out); end
fprintf(fid, 'instance,f_target,source,date\n');
for r = 1:numel(T), fprintf(fid, '%s,%.16g,%s,%s\n', T(r).instance, T(r).f_target, T(r).source, T(r).date); end
fclose(fid);
fprintf('%s: %d targets (%d certified, %d f_star, %d best_known, %d none) -> %s\n', fname, numel(T), ...
    sum(strcmp({T.source}, 'certified')), sum(strcmp({T.source}, 'f_star')), sum(strcmp({T.source}, 'best_known')), sum(strcmp({T.source}, 'none')), out);
end

function S = read_csv(file)
% header + rows -> struct of cell columns (all strings)
S = struct();
fid = fopen(file, 'r'); hdr = strsplit(strtrim(fgetl(fid)), ',');
for k = 1:numel(hdr), hdr{k} = regexprep(hdr{k}, '\W', '_'); S.(hdr{k}) = {}; end
while true
    line = fgetl(fid); if ~ischar(line), break; end
    c = strsplit(line, ',');
    for k = 1:numel(hdr)
        if k <= numel(c), S.(hdr{k}){end+1} = c{k}; else, S.(hdr{k}){end+1} = ''; end
    end
end
fclose(fid);
end
