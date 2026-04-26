function [lb_scaled, ub_scaled, info] = generate_scaled_bounds(lb_orig, ub_orig, scale_factors, varargin)
% GENERATE_SCALED_BOUNDS - Generate scaled bounds with numerical safety
%
% Syntax:
%   [lb_s, ub_s] = generate_scaled_bounds(lb, ub, s)
%   [lb_s, ub_s, info] = generate_scaled_bounds(lb, ub, s)
%   [lb_s, ub_s, info] = generate_scaled_bounds(lb, ub, s, 'AsymmetryRatio', alpha)
%
% Inputs:
%   lb_orig        - Original lower bounds [n×1]
%   ub_orig        - Original upper bounds [n×1]
%   scale_factors  - Scale factors [n×1] with s_i > 0
%
% Optional Name-Value Parameters:
%   'AsymmetryRatio' - α ∈ [0,1], fraction of scaled width above center
%                      (default: 0.5 for symmetric scaling)
%
% Outputs:
%   lb_scaled - Scaled lower bounds [n×1]
%   ub_scaled - Scaled upper bounds [n×1]
%   info      - Structure with transformation details
%
% Numerical Safety:
%   MIN_WIDTH = 1e-3 - Minimum amplitude for scaled intervals
%   If width_scaled(i) = scale_factors(i) × (ub_orig(i) - lb_orig(i)) < MIN_WIDTH,
%   then scale_factors(i) is adjusted to ensure width_scaled(i) = MIN_WIDTH.
%   This prevents division by small values in coordinate transformations.
%
% Transformation:
%   For symmetric scaling (α = 0.5, default):
%     center_orig(i) = (lb_orig(i) + ub_orig(i)) / 2
%     center_scaled(i) = center_orig(i) × scale_factors(i)
%     width_scaled(i) = scale_factors(i) × (ub_orig(i) - lb_orig(i))
%     lb_scaled(i) = center_scaled(i) - width_scaled(i) / 2
%     ub_scaled(i) = center_scaled(i) + width_scaled(i) / 2
%
%   For asymmetric scaling (α ≠ 0.5):
%     lb_scaled(i) = center_scaled(i) - α × width_scaled(i)
%     ub_scaled(i) = center_scaled(i) + (1-α) × width_scaled(i)
%
% Examples:
%   % Symmetric scaling (preserves center)
%   lb = [0; 0];
%   ub = [10; 100];
%   s = [1; 10];
%   [lb_s, ub_s] = generate_scaled_bounds(lb, ub, s);
%   % Result: [0, 10] and [0, 1000]
%
%   % Asymmetric scaling (70% above, 30% below center)
%   [lb_s, ub_s] = generate_scaled_bounds(lb, ub, s, 'AsymmetryRatio', 0.7);
%
%   % Get detailed info
%   [lb_s, ub_s, info] = generate_scaled_bounds(lb, ub, s);
%   fprintf('Contrast: %.2e\n', info.contrast);
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

%% Numerical safety threshold
MIN_WIDTH = 1e-3;  % Minimum amplitude of scaled interval

%% Input validation
narginchk(3, inf);

% Ensure column vectors FIRST
lb_orig = lb_orig(:);
ub_orig = ub_orig(:);
scale_factors = scale_factors(:);

% Parse optional parameters
p = inputParser;
addParameter(p, 'AsymmetryRatio', 0.5, @(x) isscalar(x) && x >= 0 && x <= 1);
parse(p, varargin{:});
opts = p.Results;

n = length(lb_orig);

% Validate dimensions
if length(ub_orig) ~= n || length(scale_factors) ~= n
    error('generate_scaled_bounds:DimensionMismatch', ...
          'All inputs must have same dimension. Got: lb(%d), ub(%d), s(%d)', ...
          length(lb_orig), length(ub_orig), length(scale_factors));
end

% Validate that inputs are numeric
if ~isnumeric(lb_orig) || ~isnumeric(ub_orig) || ~isnumeric(scale_factors)
    error('generate_scaled_bounds:InvalidInputType', ...
          'All inputs must be numeric vectors.');
end

% Validate bounds
if any(ub_orig <= lb_orig)
    error('generate_scaled_bounds:InvalidBounds', ...
          'Upper bounds must be strictly greater than lower bounds.');
end

% Validate scale factors
if any(scale_factors <= 0)
    error('generate_scaled_bounds:InvalidScaleFactors', ...
          'Scale factors must be > 0. Got min: %.2e', min(scale_factors));
end

%% Apply numerical safety: ensure minimum width
width_orig = ub_orig - lb_orig;
width_scaled = scale_factors .* width_orig;

% Store original scale factors
scale_factors_original = scale_factors;

% Find dimensions where scaled width is too small
too_small = width_scaled < MIN_WIDTH;

if any(too_small)
    % Adjust scale factors to ensure minimum width
    % width_scaled(i) = scale_factors(i) × width_orig(i) ≥ MIN_WIDTH
    % => scale_factors(i) = MIN_WIDTH / width_orig(i)
    scale_factors(too_small) = MIN_WIDTH ./ width_orig(too_small);
    
    % Recalculate widths
    width_scaled = scale_factors .* width_orig;
    
    % Report adjustment
    n_adjusted = sum(too_small);
    idx_adjusted = find(too_small);
    
    fprintf('⚠ Numerical safety adjustment applied:\n');
    fprintf('  MIN_WIDTH threshold: %.2e\n', MIN_WIDTH);
    fprintf('  Dimensions adjusted: %d/%d → [', n_adjusted, n);
    fprintf('%d ', idx_adjusted);
    fprintf(']\n');
    
    % Show example
    if n_adjusted > 0
        i_ex = idx_adjusted(1);
        fprintf('  Example (dim %d):\n', i_ex);
        fprintf('    Original:  s=%.2e, w_orig=%.2e → w_scaled=%.2e\n', ...
                scale_factors_original(i_ex), width_orig(i_ex), ...
                scale_factors_original(i_ex) * width_orig(i_ex));
        fprintf('    Adjusted:  s=%.2e, w_orig=%.2e → w_scaled=%.2e\n', ...
                scale_factors(i_ex), width_orig(i_ex), width_scaled(i_ex));
    end
    fprintf('\n');
end

%% Compute scaled bounds (preserve center with asymmetry)
% CRITICAL FIX: Scale the center, not just the width!
center_orig = (lb_orig + ub_orig) / 2;
center_scaled = center_orig .* scale_factors;  % Scale the center!
alpha = opts.AsymmetryRatio;

% Symmetric scaling formula (default α = 0.5):
% lb_scaled = center_scaled - α × width_scaled
% ub_scaled = center_scaled + (1-α) × width_scaled
lb_scaled = center_scaled - alpha * width_scaled;
ub_scaled = center_scaled + (1 - alpha) * width_scaled;

%% Build info struct
if nargout >= 3
    info = struct();
    
    % Basic dimensions
    info.n = n;
    info.scale_factors = scale_factors;
    info.scale_factors_original = scale_factors_original;
    info.scale_factors_adjusted = any(too_small);
    info.n_adjusted = sum(too_small);
    info.idx_adjusted = find(too_small)';
    
    % Original bounds
    info.lb_orig = lb_orig;
    info.ub_orig = ub_orig;
    info.width_orig = width_orig;
    info.center_orig = center_orig;
    
    % Scaled bounds
    info.lb_scaled = lb_scaled;
    info.ub_scaled = ub_scaled;
    info.width_scaled = width_scaled;
    info.center_scaled = center_scaled;  % Use the scaled center
    
    % Contrast (using adjusted scale factors)
    info.contrast = max(width_scaled) / min(width_scaled);
    info.contrast_original = max(scale_factors_original) / min(scale_factors_original);
    
    % Asymmetry parameter
    info.asymmetry_ratio = alpha;
    
    % Numerical safety
    info.min_width_threshold = MIN_WIDTH;
    info.min_width_actual = min(width_scaled);
    
    % Verify center preservation (scaled centers should match)
    if abs(alpha - 0.5) < 1e-10
        % For symmetric scaling, verify that:
        % center_scaled = (lb_scaled + ub_scaled) / 2
        computed_center = (lb_scaled + ub_scaled) / 2;
        center_error = norm(computed_center - center_scaled);
        info.center_preservation_error = center_error;
        
        if center_error > 1e-10
            warning('generate_scaled_bounds:CenterNotPreserved', ...
                    'Center not preserved! Error: %.2e', center_error);
        end
    end
    
    % Statistics
    info.min_scale = min(scale_factors);
    info.max_scale = max(scale_factors);
    info.mean_scale = mean(scale_factors);
    info.std_scale = std(scale_factors);
end

%% Display summary if no output requested
if nargout == 0
    display_summary(lb_orig, ub_orig, lb_scaled, ub_scaled, scale_factors, alpha, n);
end

end

%% Helper function: Display summary
function display_summary(lb_orig, ub_orig, lb_scaled, ub_scaled, scale_factors, alpha, n)
    
    fprintf('\n=== Scaled Bounds Generated ===\n');
    fprintf('Dimension: n = %d\n', n);
    fprintf('Contrast:  κ = %.2e\n', max(scale_factors)/min(scale_factors));
    
    if abs(alpha - 0.5) < 1e-10
        fprintf('Mode:      Symmetric (α = 0.5, center preserved)\n');
    else
        fprintf('Mode:      Asymmetric (α = %.2f)\n', alpha);
        fprintf('           %.0f%% above center, %.0f%% below\n', ...
                alpha*100, (1-alpha)*100);
    end
    
    width_orig = ub_orig - lb_orig;
    width_scaled = ub_scaled - lb_scaled;
    
    fprintf('\nWidth statistics:\n');
    fprintf('  Original:  [%.2e, %.2e]\n', min(width_orig), max(width_orig));
    fprintf('  Scaled:    [%.2e, %.2e]\n', min(width_scaled), max(width_scaled));
    fprintf('  Expansion: %.1fx to %.1fx\n', ...
            min(width_scaled./width_orig), max(width_scaled./width_orig));
    
    if n <= 10
        fprintf('\nDimension-wise details:\n');
        fprintf('%3s | %12s | %12s | %12s | %12s | %10s | %10s\n', ...
                'i', 'lb_orig', 'ub_orig', 'lb_scaled', 'ub_scaled', 'scale', 'width_sc');
        fprintf('%s\n', repmat('-', 1, 95));
        
        for i = 1:n
            fprintf('%3d | %12.4f | %12.4f | %12.4f | %12.4f | %10.2e | %10.2e\n', ...
                    i, lb_orig(i), ub_orig(i), lb_scaled(i), ub_scaled(i), ...
                    scale_factors(i), width_scaled(i));
        end
    else
        fprintf('\n(Showing first 5 and last 5 dimensions)\n');
        fprintf('%3s | %12s | %12s | %12s | %12s | %10s | %10s\n', ...
                'i', 'lb_orig', 'ub_orig', 'lb_scaled', 'ub_scaled', 'scale', 'width_sc');
        fprintf('%s\n', repmat('-', 1, 95));
        
        % First 5
        for i = 1:min(5, n)
            fprintf('%3d | %12.4f | %12.4f | %12.4f | %12.4f | %10.2e | %10.2e\n', ...
                    i, lb_orig(i), ub_orig(i), lb_scaled(i), ub_scaled(i), ...
                    scale_factors(i), width_scaled(i));
        end
        
        if n > 10
            fprintf('... | %12s | %12s | %12s | %12s | %10s | %10s\n', ...
                    '...', '...', '...', '...', '...', '...');
        end
        
        % Last 5
        for i = max(n-4, 6):n
            fprintf('%3d | %12.4f | %12.4f | %12.4f | %12.4f | %10.2e | %10.2e\n', ...
                    i, lb_orig(i), ub_orig(i), lb_scaled(i), ub_scaled(i), ...
                    scale_factors(i), width_scaled(i));
        end
    end
    
    fprintf('\n');
end