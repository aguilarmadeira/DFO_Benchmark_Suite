function H = holdout_instances()
%HOLDOUT_INSTANCES  The held-out part of the DFO Benchmark Suite v2 core: for validation only.
%
%   H = holdout_instances()      H.base_instances, H.is_holdout(PB), H.rule
%
% Convention adopted on 2026-11-06, before any solver was tuned on the version-2 collection:
% results on the held-out part are reported (separately from the development part), but no
% parameter of a solver and no variant of an algorithm is chosen using them.
%
% The held-out part is defined by base function and dimension, so it follows the base instance
% across the type patterns (C, D, K, DK, CD, CK, CDK) and scale strategies: an instance of the core
% is held out when PB.base_function and PB.nvar match one of the 19 base instances below (every
% third of the 63 original base instances, in their original order; 2 of those 21 are not in the
% core). Family BK (folders BK_MK16, BK_binary) has no held-out part in version 2.0.0: it contains
% seeds 1-3 only, which the convention assigns to development.
H.base_instances = {'becker_lago_2D', 'cosine_mixture_4D', 'fifteenn_local_minima_2D', ...
    'fifteenn_local_minima_8D', 'goldstein_price_2D', 'gulf_3D', 'hosaki_2D', 'mccormick_2D', ...
    'multi_gaussian_2D', 'paviani_10D', 'powell_4D', 'salomon_5D', 'schaffer2_2D', 'shekel_47_4D', ...
    'shekel_foxholes_10D', 'sixhumpcamel_2D', 'tenn_local_minima_2D', 'tenn_local_minima_8D', 'wood_4D'};
H.is_holdout = @(PB) isfield(PB, 'base_function') && ...
    any(strcmp(sprintf('%s_%dD', PB.base_function, PB.nvar), H.base_instances));
H.rule = ['Held out: the 19 base instances listed (by base_function and dimension), in every pattern ' ...
          'and strategy folder; Family BK: none in v2.0.0. Not to be used for parameter selection or algorithm development.'];
end
