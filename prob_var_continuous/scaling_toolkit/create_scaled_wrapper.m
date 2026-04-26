function fwrap = create_scaled_wrapper(fhandle, lb_orig, ub_orig, lb_work, ub_work)
%CREATE_SCALED_WRAPPER  Build a scaled-domain wrapper for a bound-constrained objective.
%
%   FWRAP = CREATE_SCALED_WRAPPER(FHANDLE, LB_ORIG, UB_ORIG, LB_WORK, UB_WORK)
%   returns a function handle FWRAP such that FWRAP(X_WORK) evaluates
%   FHANDLE at the point in the original domain corresponding to the same
%   relative position within [LB_WORK, UB_WORK]. Specifically:
%
%       t      = clip01( (x_work - lb_work) ./ (ub_work - lb_work) )
%       x_orig = lb_orig + t .* (ub_orig - lb_orig)
%       fwrap(x_work) = fhandle(x_orig)
%
%   This implements equation (B2) of the accompanying paper.
%
%   BOUNDARY HANDLING (clip-to-feasible):
%   The mapping clips infeasible inputs componentwise to [0,1] before
%   denormalization, mirroring the convention used in the 504 self-contained
%   wrappers shipped with this benchmark suite (see info.mapping field):
%
%       t = clip01( (x_work - lb_work) ./ (ub_work - lb_work) )
%
%   This guarantees that fwrap(x_work) is always finite and well-defined,
%   even for x_work slightly outside [lb_work, ub_work] (e.g., due to
%   floating-point boundary effects). For algorithms that require an
%   extreme-barrier interface — i.e., +Inf at infeasible points — this
%   should be enforced at the algorithm level (as is done internally by
%   GLODS-SI), not at the wrapper level.
%
%   Inputs:
%       fhandle  - Function handle of the original objective f : R^n -> R
%       lb_orig  - Original lower bounds [n x 1]
%       ub_orig  - Original upper bounds [n x 1]
%       lb_work  - Working (scaled) lower bounds [n x 1]
%       ub_work  - Working (scaled) upper bounds [n x 1]
%
%   Output:
%       fwrap    - Function handle, callable as fwrap(x_work).
%
%   Example:
%       % Original objective on [0,1]^2
%       f       = @(x) sum(x.^2);
%       lb_o    = [0; 0];
%       ub_o    = [1; 1];
%       % Scaled work domain (10x in second coordinate)
%       lb_w    = [0;  0];
%       ub_w    = [1; 10];
%       fwrap   = create_scaled_wrapper(f, lb_o, ub_o, lb_w, ub_w);
%       fwrap([0.5; 5])    % evaluates f at relative center, i.e. f([0.5; 0.5])
%
% Reference:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.
%
% Copyright (C) 2026 J. F. A. Madeira.
% SPDX-License-Identifier: LGPL-3.0-or-later

% Ensure column vectors.
lb_orig = lb_orig(:);
ub_orig = ub_orig(:);
lb_work = lb_work(:);
ub_work = ub_work(:);

% Validate dimensions.
n = numel(lb_orig);
if numel(ub_orig) ~= n || numel(lb_work) ~= n || numel(ub_work) ~= n
    error('create_scaled_wrapper:DimensionMismatch', ...
        'All bound vectors must have the same length.');
end

% Validate ordering.
if any(ub_orig <= lb_orig)
    error('create_scaled_wrapper:InvalidOriginalBounds', ...
        'ub_orig must be strictly greater than lb_orig.');
end
if any(ub_work <= lb_work)
    error('create_scaled_wrapper:InvalidWorkBounds', ...
        'ub_work must be strictly greater than lb_work.');
end

% Pre-compute affine map coefficients (captured by the closure below).
range_work = ub_work - lb_work;
range_orig = ub_orig - lb_orig;

% Build the wrapper.
%   x_work   ---->  t = clip01((x_work - lb_work) ./ range_work)
%   x_orig   =  lb_orig + t .* range_orig
%   fwrap(x_work) = fhandle(x_orig)
fwrap = @(x_work) fhandle( ...
    lb_orig + ...
    max(0, min(1, (x_work(:) - lb_work) ./ range_work)) .* range_orig );
end
