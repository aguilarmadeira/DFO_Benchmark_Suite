function [P, rec] = bench_problem(PB, budget, seed, save_points)
%BENCH_PROBLEM  Counted view of a suite problem, as handed to a solver.
%
%   [P, rec] = bench_problem(PB, budget, seed)
%
% P is a struct with
%   P.name, P.nvar, P.var_types, P.cardinality (Inf for continuous),
%   P.ProblemData, P.lb_C_work / P.ub_C_work (continuous bounds, if any),
%   P.func_F(xr)   counted objective, xr cell of values/labels (as PB.func_F)
%   P.func_idx(z)  counted objective, z vector of 1-based indices
%                  (continuous coordinates as values)         (as PB.func_idx)
%   P.idx2cell, P.cell2idx  conversions (not counted)
%   P.x0_idx       the reference starting point of the protocol: PB.x0_idx when
%                  the wrapper carries one (version-2 core: a Halton point), else
%                  the canonical centre (index round(0.5*(m-1))+1 per finite
%                  variable, midpoint of the bounds for continuous variables)
%   P.x0           the same point as a cell
%   P.budget, P.seed
% rec is the BenchRecorder that holds the trace (see run_benchmark).
%
% The solver must use ONLY P.func_F or P.func_idx to evaluate; everything else
% about the problem (f_star, x_star, ...) is deliberately not exposed.
if nargin < 2 || isempty(budget), budget = 1000 * PB.nvar; end
if nargin < 3 || isempty(seed), seed = 1; end
if nargin < 4 || isempty(save_points), save_points = false; end
rec = BenchRecorder(budget, save_points);
P = struct();
P.name        = PB.name;
P.nvar        = PB.nvar;
P.var_types   = PB.var_types;
P.ProblemData = PB.ProblemData;
if isfield(PB, 'cardinality'), P.cardinality = PB.cardinality; else, P.cardinality = local_cardinality(PB); end
if isfield(PB, 'lb_C_work'), P.lb_C_work = PB.lb_C_work; P.ub_C_work = PB.ub_C_work; end
% both interfaces are counted, and the DISTINCT points are counted on a common
% key (the index representation), so that a point evaluated through func_F and
% through func_idx is the same point
if isfield(PB, 'cell2idx')
    P.func_F   = @(xr) rec.record(PB.func_F(xr), xr, point_key(PB.cell2idx(xr)), PB.cell2idx(xr));
else
    P.func_F   = @(xr) rec.record(PB.func_F(xr), xr);
end
if isfield(PB, 'func_idx')
    P.func_idx = @(z) rec.record(PB.func_idx(z), z, point_key(z), z);
    P.idx2cell = PB.idx2cell;
    P.cell2idx = PB.cell2idx;
else
    P.func_idx = @(z) error('bench_problem: this wrapper (version 1) has no index interface; use P.func_F.');
end
n = PB.nvar; z0 = zeros(1, n); x0 = cell(1, n);
if isfield(PB, 'x0_idx')                      % version-2 core: the Halton starting point of the protocol
    z0 = PB.x0_idx(:)';
    for i = 1:n
        switch PB.var_types{i}
            case 'C', x0{i} = z0(i);
            case 'D', x0{i} = PB.ProblemData{i}{1}(z0(i));
            case 'K', x0{i} = PB.ProblemData{i}{z0(i)};
        end
    end
else                                          % earlier folders: the canonical centre
    for i = 1:n
        switch PB.var_types{i}
            case 'C'
                z0(i) = 0.5 * (PB.lb_C_work(i) + PB.ub_C_work(i)); x0{i} = z0(i);
            case 'D'
                m = P.cardinality(i); z0(i) = round(0.5*(m-1)) + 1; x0{i} = PB.ProblemData{i}{1}(z0(i));
            case 'K'
                m = P.cardinality(i); z0(i) = round(0.5*(m-1)) + 1; x0{i} = PB.ProblemData{i}{z0(i)};
        end
    end
end
P.x0_idx = z0; P.x0 = x0;
P.budget = budget; P.seed = seed;
end

function k = point_key(z)
k = sprintf('%.17g,', z);
end

function c = local_cardinality(PB)
n = PB.nvar; c = zeros(1, n);
for i = 1:n
    pd = PB.ProblemData{i};
    if strcmp(PB.var_types{i}, 'C'), c(i) = Inf;
    elseif iscell(pd) && isnumeric(pd{1}), c(i) = numel(pd{1});
    else, c(i) = numel(pd);
    end
end
end
