function x = halton_point(idx, d)
%HALTON_POINT  Points idx (1-based, vector) of the d-dimensional Halton sequence, bases = first d primes.
%   x = halton_point(idx, d)   -> numel(idx) x d, entries in (0,1)
% Exact radical inverse evaluated with the same floating-point operations in MATLAB and Octave
% (sums of i_k / b^k accumulated from the most significant digit downwards would differ by
% rounding order; this implementation accumulates in a fixed order, so both give identical
% doubles — verified on the published files). Stateless: no random generator,
% no seed. The suite's deterministic construction principle uses it for starting points and
% for every choice that needs well-distributed values (categorical permutations through a
% Fisher-Yates shuffle driven by consecutive coordinates).
P = primes(max(30, 20 * d)); P = P(1:d); idx = idx(:); x = zeros(numel(idx), d);
for j = 1:d
    b = P(j);
    for r = 1:numel(idx)
        i = idx(r); f = 1 / b; v = 0;
        while i > 0
            v = v + f * mod(i, b); i = floor(i / b); f = f / b;
        end
        x(r, j) = v;
    end
end
end
