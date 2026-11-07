function PB = attach_reference_optimum(PB, folder)
%ATTACH_REFERENCE_OPTIMUM  Add the optimal solution of a version-2 core instance to PB (for validation and
% path statistics only; solvers never receive it). The core wrappers do not carry it: it is read from the
% folder's reference_optima.csv (instance, base_function, n, f_star, justification, n_min, x_star_idx, x_star_orig)
if isfield(PB, 'x_star_all_idx'), return; end
file = fullfile(folder, 'reference_optima.csv');
if exist(file, 'file') ~= 2, return; end
fid = fopen(file, 'r'); z = [];
while true
    line = fgetl(fid); if ~ischar(line), break; end
    c = strsplit(strtrim(line), ','); if numel(c) < 7, continue; end
    if strcmp(c{1}, PB.name), z = str2double(strsplit(c{7}, ';')); break; end
end
fclose(fid);
if isempty(z), return; end
PB.x_star_idx = z(:)'; PB.x_star_all_idx = {PB.x_star_idx};
PB.x_star = PB.idx2cell(PB.x_star_idx); PB.x_star_all = {PB.x_star};
if ~isfield(PB, 'f_star_documented'), PB.f_star_documented = PB.f_star; end
if ~isfield(PB, 'domain_size'), PB.domain_size = NaN; end
end
