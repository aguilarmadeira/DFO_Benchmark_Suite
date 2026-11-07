function perm = halton_permutation(m, idx)
%HALTON_PERMUTATION  Fisher-Yates shuffle of 1:m driven by the Halton point of index idx in
% dimension m-1 (one coordinate per step). Deterministic and portable; reproducible from (m, idx).
% Not a uniformly distributed permutation in the probabilistic sense (the Halton coordinates are
% not independent), but free of a systematic pattern and identical in MATLAB and Octave.
perm = (1:m)';
if m < 2, return; end
u = halton_point(idx, m - 1);
for i = m:-1:2
    j = 1 + floor(u(m - i + 1) * i); if j > i, j = i; end
    t = perm(i); perm(i) = perm(j); perm(j) = t;
end
end
