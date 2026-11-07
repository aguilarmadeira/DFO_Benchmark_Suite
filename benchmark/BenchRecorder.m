classdef BenchRecorder < handle
%BENCHRECORDER  Evaluation counter and trace for one benchmark run.
%
% Created by bench_problem; solvers never see it directly (they receive the
% counted handles P.func_F / P.func_idx). It stores the number of
% evaluations, the best value, the first point that attained it, and the
% trace of improvements (evaluation index, best value so far). When the
% budget is exhausted every further evaluation raises the error
% 'BENCH:budget', which run_benchmark catches: a solver that ignores the
% budget is stopped without penalty, a solver that respects it never sees it.
% A non-finite value (NaN, +-Inf: e.g. meyer_roth has a singular plane inside
% the box) counts as one evaluation and is returned to the solver as +Inf, so
% every solver receives the same information about an invalid point.
    properties
        budget   = Inf;
        nevals   = 0;
        fbest    = Inf;
        xbest    = [];       % cell (func_F) or index vector (func_idx) of the first best point
        trace    = zeros(0, 2);   % [neval fbest] at every strict improvement
        n_nonfinite = 0;     % evaluations that returned NaN/Inf
        keys = {};           % key (index vector as text) of every evaluation, in order
        nk = 0;              % number of keys stored
        save_points = false; % when true, every evaluated point (index representation) and its value are kept
        points = zeros(0, 0);% [z f] per evaluation, in order (only with save_points)
    end
    properties (Dependent)
        n_distinct           % number of distinct points evaluated (a solver without cache may re-evaluate x0)
    end
    methods
        function v = get.n_distinct(obj)
            v = numel(unique(obj.keys(1:obj.nk)));       % computed on demand: O(N log N), no map (fast in Octave too)
        end
        function obj = BenchRecorder(budget, save_points)
            if nargin >= 1 && ~isempty(budget), obj.budget = budget; end
            if nargin >= 2 && ~isempty(save_points), obj.save_points = save_points; end
        end
        function f = record(obj, f, x, key, z)
            if obj.nevals >= obj.budget
                error('BENCH:budget', 'evaluation budget (%d) exhausted', obj.budget);
            end
            obj.nevals = obj.nevals + 1;
            if nargin >= 4 && ~isempty(key)
                if obj.nk >= numel(obj.keys), obj.keys{max(1024, 2 * numel(obj.keys))} = ''; end   % grow geometrically
                obj.nk = obj.nk + 1; obj.keys{obj.nk} = key;
            end
            if ~isfinite(f), obj.n_nonfinite = obj.n_nonfinite + 1; f = Inf; end
            if obj.save_points && nargin >= 5 && ~isempty(z)
                if isempty(obj.points), obj.points = NaN(min(obj.budget, 1e6), numel(z) + 1); end
                if obj.nevals <= size(obj.points, 1), obj.points(obj.nevals, :) = [z(:)', f]; end
            end
            if f < obj.fbest
                obj.fbest = f;
                obj.xbest = x;
                obj.trace(end+1, :) = [obj.nevals, f];
            end
        end
        function reset(obj)
            obj.nevals = 0; obj.fbest = Inf; obj.xbest = []; obj.trace = zeros(0,2); obj.n_nonfinite = 0; obj.keys = {}; obj.nk = 0; obj.points = zeros(0,0);
        end
    end
end
