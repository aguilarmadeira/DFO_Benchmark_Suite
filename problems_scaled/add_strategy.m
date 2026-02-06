function add_strategy(strategy_name)
%ADD_STRATEGY Add one strategy folder under ./problemsNew to the MATLAB path.
% Usage: add_strategy('progressive')

root = fileparts(mfilename('fullpath'));
s = lower(strategy_name);
switch s
    case 'baseline'
        addpath(fullfile(root, 'baseline_1'));
    case 'progressive'
        addpath(fullfile(root, 'progressive_1e6'));
    case 'extreme'
        addpath(fullfile(root, 'extreme_1e8'));
    case 'sobol_oscillatory'
        addpath(fullfile(root, 'sobol_oscillatory_1e6'));
    case 'sobol_digit_oscillatory'
        addpath(fullfile(root, 'sobol_digit_oscillatory_1e8'));
    case 'spatial_thermal'
        addpath(fullfile(root, 'spatial_thermal_9e4'));
    case 'halton_oscillatory'
        addpath(fullfile(root, 'halton_oscillatory_1e6'));
    otherwise
        error('Unknown strategy: %s', strategy_name);
end
end
