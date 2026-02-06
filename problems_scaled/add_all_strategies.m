function add_all_strategies()
%ADD_ALL_STRATEGIES Add all strategy folders under ./problemsNew to the MATLAB path.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'baseline_1'));
addpath(fullfile(root, 'progressive_1e6'));
addpath(fullfile(root, 'extreme_1e8'));
addpath(fullfile(root, 'sobol_oscillatory_1e6'));
addpath(fullfile(root, 'sobol_digit_oscillatory_1e8'));
addpath(fullfile(root, 'spatial_thermal_9e4'));
addpath(fullfile(root, 'halton_oscillatory_1e6'));
end
