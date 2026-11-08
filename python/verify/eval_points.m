% eval_points.m — evaluates the MATLAB wrappers of the DFO Benchmark Suite v2.1.0 at the points of points.csv
% (x0, the 9 alternative starting points, x* and 20 random points per instance), to compare with the Python version.
% Location: python/verify of the suite; the MATLAB wrappers are read from core/ at the root of the suite.
% Output: values_matlab.csv in this folder; compare with  python compare_values.py
here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
fid = fopen(fullfile(here, 'points.csv'), 'r'); out = fopen(fullfile(here, 'values_matlab.csv'), 'w');
cur = ''; PB = []; t0 = tic; k = 0;
while true
    l = fgetl(fid); if ~ischar(l), break; end
    c = strsplit(l, ',');
    if ~strcmp(cur, [c{1} '/' c{2}])
        if ~isempty(cur), rmpath(fullfile(root, 'core', fo)); end
        fo = c{1}; addpath(fullfile(root, 'core', fo));
        try, [~, PB] = evalc(['feval(''' c{2} ''')']); catch, evalc(['PB = feval(''' c{2} ''');']); end   % MATLAB / Octave
        cur = [c{1} '/' c{2}];
    end
    z = str2double(strsplit(c{4}, ';'));
    fprintf(out, '%s,%s,%s,%.17g\n', c{1}, c{2}, c{3}, PB.func_idx(z));
    k = k + 1; if mod(k, 5000) == 0, fprintf('%d points, %.0f s\n', k, toc(t0)); end
end
fclose(fid); fclose(out); rmpath(fullfile(root, 'core', fo));
fprintf('done: %d points in %.0f s -> values_matlab.csv\n', k, toc(t0));
