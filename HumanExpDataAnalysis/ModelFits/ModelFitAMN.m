close all
clear all
restoredefaultpath

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsModelAMN')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsAMN')

% expData            = load('Data\CORNFB01.mat'); % Yichao
% expData            = load('./Data/COR33.mat'); % Akash
% expData            = load('./Data/COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming
expData            = load('./Data/CORNFB02.mat');% Jonathan

fltData       = expData.dat( expData.dat.session > 2 , :); 
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors
initCond      = getInitialConditions(formattedData);

%%
optParams.nStarts = 30;
result = Optimize(formattedData, initCond, [], optParams);

%%
save("./Data/jonathan_fit_model_AMN.mat", "result");

%%
% load('./Data/jonathan_fit_model_AMN.mat');
errBins   = -90:3:90; 

[~, idx] = min(result.f);

n_uncertainty_levels = 6;

opt_param_sigma_ext         = result.x(idx, 1:n_uncertainty_levels);
opt_param_ext_noise_scale   = result.x(idx, n_uncertainty_levels+1:2*n_uncertainty_levels);
opt_param_scale             = result.x(idx, 2*n_uncertainty_levels + 1);
opt_param_sigma_meta        = result.x(idx, 2*n_uncertainty_levels + 2);
opt_param_Cc                = result.x(idx, 2*n_uncertainty_levels + 3);
opt_param_guessrate         = result.x(idx, 2*n_uncertainty_levels + 4);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("sigma_ext: Fit: %.4f C: %.3f, S: %.3f, D: %.3f\n", ...
        opt_param_sigma_ext(i), ...
        formattedData.uncertaintyVals(i, 1), ...
        formattedData.uncertaintyVals(i, 2), ...
        formattedData.uncertaintyVals(i, 3))
end

for i =1:n_uncertainty_levels
    fprintf("sigma_ext scale Fit: %.4f \n", opt_param_ext_noise_scale(i))
end

fprintf("Scale Fit: %.4f \n", opt_param_scale)
fprintf("Meta Fit: %.4f \n", opt_param_sigma_meta)
fprintf("Cc Fit: %.4f \n", opt_param_Cc)
fprintf("GR Fit: %.4f \n", opt_param_guessrate)

for i =1:n_uncertainty_levels
    fprintf("\n %.4f C: %.3f, S: %.3f, D: %.3f \n", ...
        opt_param_sigma_ext(i)*opt_param_ext_noise_scale(i), ...
        formattedData.uncertaintyVals(i, 1), ...
        formattedData.uncertaintyVals(i, 2), ...
        formattedData.uncertaintyVals(i, 3));
end

% full model - might not be nice for human subjects - might need more
% customization

%% TODO: plot
modelParams = result.x(idx, :);
plotFitResultAMN(formattedData, modelParams, initCond, errBins)

