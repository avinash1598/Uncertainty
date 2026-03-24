close all
clear all
restoredefaultpath

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsModelAMN_v2')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsAMN_v2')

% expData            = load('Data\CORNFB01.mat'); % Yichao
expData            = load('./Data/COR33.mat'); % Akash
% expData            = load('./Data/COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming
% expData            = load('./Data/CORNFB02.mat');% Jonathan

fltData       = expData.dat( expData.dat.session > 0 , :); 
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors
initCond      = getInitialConditions(formattedData);

%%
optParams.nStarts = 30;
result = Optimize(formattedData, initCond, [], optParams);

%%
save("./Data/akash_fit_model_AMN_v2.mat", "result");

%%
% load('./Data/akash_fit_model_AMN_v2.mat');
errBins   = -90:3:90; 

[~, idx] = min(result.f);

n_uncertainty_levels = 6;

opt_param_sigma_ext         = result.x(idx, 1:n_uncertainty_levels);
opt_param_scale             = result.x(idx, n_uncertainty_levels + 1);
opt_param_sigma_meta        = result.x(idx, n_uncertainty_levels + 2);
opt_param_Cc                = result.x(idx, n_uncertainty_levels + 3);
opt_param_guessrate         = result.x(idx, n_uncertainty_levels + 4);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("sigma_ext: Fit: %.4f \n", opt_param_sigma_ext(i))
end

fprintf("Scale Fit: %.4f \n", opt_param_scale)
fprintf("Meta Fit: %.4f \n", opt_param_sigma_meta)
fprintf("Cc Fit: %.4f \n", opt_param_Cc)
fprintf("GR Fit: %.4f \n", opt_param_guessrate)

% full model - might not be nice for human subjects - might need more
% customization

%% TODO: plot
modelParams = result.x(idx, :);
plotFitResultAMN_v2(formattedData, modelParams, initCond, errBins)

