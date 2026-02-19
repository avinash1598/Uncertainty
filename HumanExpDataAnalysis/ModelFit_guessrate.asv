close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')

% expData            = load('Data\CORNFB01.mat'); % Yichao
expData            = load('./Data/COR33.mat'); % Akash
% expData            = load('Data\COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming

% Format data
errBins       = -90:1:90;
fltData       = expData.dat( expData.dat.session > 0 , :); 
f.dat         = fltData; 
formattedData = formatExpData(f, false, false);  % Keep sortByMAD to false (it is set to false in NLL script)

%% Optimization
result = Optimize(formattedData, errBins, "cov");

%%
[~, idx] = min(result.f);

n_uncertainty_levels     = formattedData.n_uncertainty_levels;
opt_param_sigma_s        = result.x(idx, 1:n_uncertainty_levels);
opt_param_scale          = result.x(idx, n_uncertainty_levels + 1);
opt_param_sigma_meta     = result.x(idx, n_uncertainty_levels + 2);
opt_param_Cc             = result.x(idx, n_uncertainty_levels + 3);
opt_param_guessrate      = result.x(idx, n_uncertainty_levels + 4);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("Fit: %.4f \n", opt_param_sigma_s(i))
end

% fprintf("Shape: %.4f \n", opt_param_shape)
fprintf("Scale: %.4f \n", opt_param_scale)
fprintf("sigma_meta: %.4f \n", opt_param_sigma_meta)
fprintf("Cc: %.4f \n", opt_param_Cc)
fprintf("Guess rate: %.4f \n", opt_param_guessrate)

%% Plot result 
% Cov
paramCovModel = result.x(idx, :);
plotFitResult_guessrate(formattedData, paramCovModel, "cov");



%% Optimization
result = Optimize(formattedData, errBins, "ind");

%%
[~, idx] = min(result.f);

n_uncertainty_levels     = formattedData.n_uncertainty_levels;
opt_param_sigma_s        = result.x(idx, 1:n_uncertainty_levels);
opt_param_shape          = result.x(idx, n_uncertainty_levels + 1);
opt_param_scale          = result.x(idx, n_uncertainty_levels + 2);
opt_param_sigma_meta     = result.x(idx, n_uncertainty_levels + 3);
opt_param_Cc             = result.x(idx, n_uncertainty_levels + 4);
opt_param_guessrate      = result.x(idx, n_uncertainty_levels + 5);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("Fit: %.4f \n", opt_param_sigma_s(i))
end

fprintf("Shape: %.4f \n", opt_param_shape)
fprintf("Scale: %.4f \n", opt_param_scale)
fprintf("sigma_meta: %.4f \n", opt_param_sigma_meta)
fprintf("Cc: %.4f \n", opt_param_Cc)
fprintf("Guess rate: %.4f \n", opt_param_guessrate)

%% Plot result 
% Ind
paramIndModel = result.x(idx, :);
plotFitResult_guessrate(formattedData, paramIndModel, "ind");

