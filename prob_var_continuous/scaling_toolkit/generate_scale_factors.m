function scale_factors = generate_scale_factors(n, strategy, varargin)
% GENERATE_SCALE_FACTORS - Generate scale factors for multi-scale testing
%
% Syntax:
%   s = generate_scale_factors(n, strategy)
%   s = generate_scale_factors(n, strategy, 'MaxContrast', kappa)
%   s = generate_scale_factors(n, strategy, 'OscillationMode', mode)
%
% Inputs:
%   n        - Problem dimension (positive integer)
%   strategy - Scaling strategy (string):
%              'baseline'              : κ = 1 (no scaling)
%              'moderate'              : κ = 10^4 (reduced cyclic pattern)
%              'progressive'           : κ = 10^6 (realistic regime)
%              'extreme'               : κ = 10^8 (stress test)
%              'sobol_oscillatory'     : κ = 10^6 (continuous QMC)
%              'sobol_digit_oscillatory': κ = 10^8 (binary QMC, digit-based)
%              'spatial_thermal'       : κ ≈ 9×10^4 (multiphysics)
%              'halton_oscillatory'    : κ = 10^6 (alternative QMC)
%
% Optional Name-Value Parameters:
%   'MaxContrast'      - Maximum contrast κ (default: depends on strategy)
%   'OscillationMode'  - For oscillatory strategies: 'alternate' or 'bipolar'
%                        Default: 'alternate'
%   'SobolSkip'        - Points to skip in Sobol sequence (default: 100)
%   'HaltonBase'       - Base for Halton sequence (default: 2)
%
% Outputs:
%   scale_factors - Column vector [n×1] with s_i > 0
%
% Examples:
%   % Baseline (no scaling)
%   s = generate_scale_factors(5, 'baseline');
%   
%   % Progressive scaling (realistic)
%   s = generate_scale_factors(10, 'progressive');
%   
%   % Extreme scaling (stress test)
%   s = generate_scale_factors(8, 'extreme');
%   
%   % Sobol oscillatory with alternate mode
%   s = generate_scale_factors(10, 'sobol_oscillatory');
%   
%   % Digit-driven oscillation
%   s = generate_scale_factors(10, 'sobol_digit_oscillatory');
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

%% Input validation
narginchk(2, inf);

if ~isscalar(n) || n <= 0 || floor(n) ~= n
    error('generate_scale_factors:InvalidDimension', ...
          'Dimension n must be a positive integer.');
end

valid_strategies = {'baseline', 'moderate', 'progressive', 'extreme', ...
                    'sobol_oscillatory', 'sobol_digit_oscillatory', ...
                    'spatial_thermal', 'halton_oscillatory'};

if ~ismember(lower(strategy), valid_strategies)
    error('generate_scale_factors:InvalidStrategy', ...
          'Unknown strategy: %s\nValid strategies: %s', ...
          strategy, strjoin(valid_strategies, ', '));
end

%% Parse optional parameters
p = inputParser;
addParameter(p, 'MaxContrast', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
addParameter(p, 'OscillationMode', 'alternate', ...
    @(x) ismember(lower(x), {'alternate', 'bipolar'}));
addParameter(p, 'SobolSkip', 100, @(x) isscalar(x) && x >= 0);
addParameter(p, 'HaltonBase', 2, @(x) isscalar(x) && x >= 2);
parse(p, varargin{:});
opts = p.Results;

%% Set default MaxContrast based on strategy
if isempty(opts.MaxContrast)
    switch lower(strategy)
        case 'baseline'
            kappa = 1;
        case 'moderate'
            kappa = 1e4;
        case 'progressive'
            kappa = 1e6;
        case 'extreme'
            kappa = 1e8;
        case 'spatial_thermal'
            kappa = 300^2;  % 90,000
        case 'sobol_oscillatory'
            kappa = 1e6;    % implicit kappa from formula s_i = 10^(6q-3)
        case 'halton_oscillatory'
            kappa = 1e6;    % same family as sobol_oscillatory
        case 'sobol_digit_oscillatory'
            kappa = 1e8;    % LOW/HIGH partition with sqrt(kappa) split
        otherwise
            kappa = 1e4;
    end
else
    kappa = opts.MaxContrast;
end

%% Generate scale factors based on strategy
switch lower(strategy)
    case 'baseline'
        scale_factors = generate_baseline(n);
        
    case 'moderate'
        scale_factors = generate_moderate(n, kappa);
        
    case 'progressive'
        scale_factors = generate_progressive(n, kappa);
        
    case 'extreme'
        scale_factors = generate_extreme(n, kappa);
        
    case 'sobol_oscillatory'
        scale_factors = generate_sobol_oscillatory(n, kappa, opts);
        
    case 'sobol_digit_oscillatory'
        scale_factors = generate_sobol_digit_oscillatory(n, kappa, opts.SobolSkip);
        
    case 'spatial_thermal'
        scale_factors = generate_spatial_thermal(n);
        
    case 'halton_oscillatory'
        scale_factors = generate_halton_oscillatory(n, kappa, opts);
        
    otherwise
        error('generate_scale_factors:UnimplementedStrategy', ...
              'Strategy %s not implemented.', strategy);
end

%% Validate output
assert(all(scale_factors > 0), 'All scale factors must be > 0');
assert(length(scale_factors) == n, 'Output dimension mismatch');
assert(isnumeric(scale_factors) && isreal(scale_factors), ...
       'Scale factors must be real numeric values');

end

%% ========================================================================
%% STRATEGY IMPLEMENTATIONS
%% ========================================================================

function s = generate_baseline(n)
% BASELINE: No scaling (κ=1) - Control condition
% Validates that normalization adds no overhead when not needed
s = ones(n, 1);
end

function s = generate_moderate(n, kappa)
% MODERATE: Reduced cyclic geometric pattern (κ=10⁴ default)
% Cyclic exponent assignment over a reduced set
% Paper R2 (Appendix B): e_i ∈ {0, 1, -1} (cyclic), s_i = 10^(η e_i)
% with η chosen so that max(s)/min(s) = κ.

exponent_sequence = [0, 1, -1];
n_exp = length(exponent_sequence);

% Assign exponents cyclically
exponents = zeros(n, 1);
for i = 1:n
    idx = mod(i - 1, n_exp) + 1;
    exponents(i) = exponent_sequence(idx);
end

% Choose η so that max(s)/min(s) = κ:
%   max e_i = +1, min e_i = -1, so range of e_i is 2
%   κ = 10^(η * 2)  =>  η = log10(κ) / 2
eta = log10(kappa) / 2;

s = 10.^(eta * exponents);
end

function s = generate_progressive(n, kappa)
% PROGRESSIVE: Realistic regime (κ=10⁶)
% Cyclic exponent assignment creates moderate disparities
% Representative of engineering problems with mixed physical units

exponent_sequence = [0, 1, 2, 3, -1, -2, -3];
n_exp = length(exponent_sequence);

% Assign exponents cyclically
exponents = zeros(n, 1);
for i = 1:n
    idx = mod(i - 1, n_exp) + 1;
    exponents(i) = exponent_sequence(idx);
end

% Scale to achieve desired κ
max_exp_target = 0.5 * log10(kappa);
current_max_exp = max(abs(exponent_sequence));
scaling = max_exp_target / current_max_exp;

s = 10.^(exponents * scaling);
end

function s = generate_extreme(n, kappa)
% EXTREME: Bipolar stress test
% First half: scale = 1 (baseline)
% Second half: scale = κ (maximum contrast)
% No intermediate values - worst case heterogeneity

half_n = floor(n / 2);
s = zeros(n, 1);

% First half: baseline scale
s(1:half_n) = 1;

% Second half: maximum scale
s(half_n+1:n) = kappa;
end


function s = generate_spatial_thermal(n)
% SPATIAL_THERMAL: Application-motivated (κ ≈ 90,000)
% First half: spatial coordinates (scale=1)
% Second half: thermal expansion coefficients (scale=300)
% Motivated by thermomechanical optimization

half_n = floor(n / 2);
s_spatial = ones(half_n, 1);
s_thermal = 300 * ones(n - half_n, 1);
s = [s_spatial; s_thermal];
end

function s = generate_sobol_oscillatory(n, kappa, opts)
% SOBOL_OSCILLATORY: Deterministic oscillation via Sobol sequence
% Tests scale-invariance under structured heterogeneity

% Generate Sobol sequence (deterministic)
p = sobolset(1);

sobol_vals = net(p, n + opts.SobolSkip);
sobol_vals = sobol_vals(opts.SobolSkip + 1 : opts.SobolSkip + n);

% Apply oscillation mode
if strcmpi(opts.OscillationMode, 'alternate')
    s = apply_alternate_oscillation(sobol_vals, kappa);
else
    s = apply_bipolar_oscillation(sobol_vals, kappa);
end
end

function s = generate_sobol_digit_oscillatory(n, kappa, skip_points)
% SOBOL_DIGIT_OSCILLATORY: Unpredictable but deterministic
% Digit-sum parity determines LOW/HIGH assignment
% Tests robustness to non-structured deterministic patterns

% Generate Sobol sequence
p = sobolset(1);
p = scramble(p, 'MatousekAffineOwen');

sobol_vals = net(p, n + skip_points);
sobol_vals = sobol_vals(skip_points + 1 : skip_points + n);

% Apply digit-sum parity rule
s = apply_sum_parity_rule(sobol_vals, kappa);
end

function s = generate_halton_oscillatory(n, kappa, opts)
% HALTON_OSCILLATORY: Alternative sequence validation
% Confirms results are not Sobol-specific

% Generate Halton sequence
halton_vals = generate_halton_sequence(n, opts.HaltonBase);

% Apply oscillation mode
if strcmpi(opts.OscillationMode, 'alternate')
    s = apply_alternate_oscillation(halton_vals, kappa);
else
    s = apply_bipolar_oscillation(halton_vals, kappa);
end
end

%% ========================================================================
%% OSCILLATION PATTERNS
%% ========================================================================

function s = apply_alternate_oscillation(quasi_vals, kappa)
% ALTERNATE: Strict L-H-L-H pattern
% Even indices → LOW [1, √κ], Odd indices → HIGH [√κ, κ]

n = length(quasi_vals);
s = zeros(n, 1);

for i = 1:n
    if mod(i, 2) == 0  % Even → LOW
        s(i) = 1 + (sqrt(kappa) - 1) * quasi_vals(i);
    else  % Odd → HIGH
        s(i) = sqrt(kappa) + (kappa - sqrt(kappa)) * quasi_vals(i);
    end
end
end

function s = apply_bipolar_oscillation(quasi_vals, kappa)
% BIPOLAR: First half LOW, second half HIGH

n = length(quasi_vals);
half_n = floor(n / 2);
s = zeros(n, 1);

% First half: LOW
for i = 1:half_n
    s(i) = 1 + (sqrt(kappa) - 1) * quasi_vals(i);
end

% Second half: HIGH
for i = half_n+1:n
    s(i) = sqrt(kappa) + (kappa - sqrt(kappa)) * quasi_vals(i);
end
end

function s = apply_sum_parity_rule(sobol_vals, kappa)
% DIGIT-SUM PARITY: Unpredictable deterministic oscillation
% Extract first 2 decimal digits, compute sum, use parity

n = length(sobol_vals);
s = zeros(n, 1);

for i = 1:n
    val = sobol_vals(i);
    
    % Extract first 2 decimal digits
    two_digits = floor(val * 100);
    digit1 = floor(two_digits / 10);
    digit2 = mod(two_digits, 10);
    digit_sum = digit1 + digit2;
    
    % Parity determines range
    if mod(digit_sum, 2) == 0  % EVEN → LOW
        s(i) = 1 + (sqrt(kappa) - 1) * sobol_vals(i);
    else  % ODD → HIGH
        s(i) = sqrt(kappa) + (kappa - sqrt(kappa)) * sobol_vals(i);
    end
end
end

%% ========================================================================
%% UTILITY FUNCTIONS
%% ========================================================================

function halton_seq = generate_halton_sequence(n, base)
% Generate Halton sequence with given base

halton_seq = zeros(n, 1);
for i = 1:n
    halton_seq(i) = van_der_corput(i, base);
end
end

function val = van_der_corput(index, base)
% Van der Corput sequence: base-b digit reversal

val = 0;
f = 1 / base;
i = index;

while i > 0
    val = val + f * mod(i, base);
    i = floor(i / base);
    f = f / base;
end
end