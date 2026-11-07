function solver = solver_random_search()
%SOLVER_RANDOM_SEARCH  Baseline for run_benchmark: quasi-random sampling of the domain.
% Points of the Halton sequence in dimension n, from index 2^20 * seed onwards (indices uniform on
% 1..m_i after scaling; continuous coordinates in the bounds); x0 has already been evaluated by
% run_benchmark. Deterministic given P.seed, identical in MATLAB and Octave (no RNG).
solver = struct('name', 'random_search', 'fun', @run, 'deterministic', false);
end

function run(P)
n = P.nvar;
isC = isinf(P.cardinality(:)');
start = 1048576 * P.seed;
while true
    U = halton_point(start + (1:512), n);       % one block of points per draw
    start = start + 512;
    for t = 1:size(U, 1)
        z = min(1 + floor(U(t, :) .* P.cardinality(:)'), P.cardinality(:)');
        if any(isC), z(isC) = P.lb_C_work(isC)' + U(t, isC) .* (P.ub_C_work(isC)' - P.lb_C_work(isC)'); end
        P.func_idx(z);                           % stops with BENCH:budget when the budget is spent
    end
end
end
