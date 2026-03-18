% close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis/Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')

% expData            = load('Data\CORNFB01.mat'); % Yichao
expData            = load('./Data/COR33.mat'); % Akash
% expData            = load('./Data/COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming

fltData       = expData.dat( expData.dat.session > 0 , :); 
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

%%
% Don't bin - this is not right
% errBins   = -90:2:90; % this is dx which might affect fitting. This value should be optimal. not too fine. not too coarse.

optParams.nStarts = 30;
optParams.hyperParamC1 = 0; % 10 or 100? Use 10 maybe to avoid overfitting
optParams.hyperParamC2 = 0;
optParams.randomGuessModel = true;

result = Optimize(formattedData, "cov", [], optParams, "full");

%%
save("akash_full_model_cov_fit_method_wraped_gauss.mat", "result");

%%
% load('akash_full_model_ind_fit_method_2.mat');
% load('akash_full_model_cov_fit_method_2.mat');
% load('tien_full_model_cov_fit_method_2.mat');
errBins   = -90:3:90; 

[~, idx] = min(result.f);

n_uncertainty_levels = 6;

opt_param_sigma_s         = result.x(idx, 1:n_uncertainty_levels);
% opt_param_shape           = result.x(idx ,n_uncertainty_levels + 1-0);
opt_param_scale           = result.x(idx ,n_uncertainty_levels + 2-1);
opt_param_sigma_meta      = result.x(idx, n_uncertainty_levels + 3-1);
opt_param_Cc              = result.x(idx, n_uncertainty_levels + 4-1);
opt_param_guessrate       = result.x(idx, n_uncertainty_levels + 5-1);
opt_param_sigma_ori_scale = result.x(idx, n_uncertainty_levels + 6-1);
opt_param_bias            = result.x(idx, n_uncertainty_levels + 7-1);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("Fit: %.4f \n", opt_param_sigma_s(i))
end

% fprintf("Shape Fit: %.4f \n", opt_param_shape)
fprintf("Scale Fit: %.4f \n", opt_param_scale)
fprintf("Meta Fit: %.4f \n", opt_param_sigma_meta)
fprintf("Cc Fit: %.4f \n", opt_param_Cc)
fprintf("GR Fit: %.4f \n", opt_param_guessrate)
fprintf("Ori scale Fit: %.4f \n", opt_param_sigma_ori_scale)
fprintf("Bias Fit: %.4f \n", opt_param_bias)

% full model - might not be nice for human subjects - might need more
% customization

%% TODO: plot
modelParams = result.x(idx, :);
% modelParams(9) = 6.8;
plotFitResult_guessrate(formattedData, modelParams, "cov", errBins, true)
