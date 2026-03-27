% TODO: do ti for all subjects

close all
clear all
restoredefaultpath

% get(groot,"default")
set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1, ...
    'defaultAxesBox','off', ...
    'defaultAxesTickDir','out', ...
    'defaultLineLineWidth',1.5, ...
    'defaultAxesLabelFontSizeMultiplier',1.1, ...
    'defaultAxesTitleFontWeight','normal', ...
    'defaultLegendBox','off');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsModel_v4/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScripts_v4/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

expData            = load('./Data/COR33.mat');     % Akash
% expData            = load('./Data/COR31.mat');   % Tien
% expData            = load('./Data/COR32.mat');   % Jiaming
% expData            = load('./Data/CORNFB02.mat');% Jonathan
% expData            = load('Data\CORNFB01.mat');  % Yichao
% For subject 6 => David: make sure to change the orientation dependent
% error

fltSessionIdx = [0, 0, 0, 0, 2, 1];
fltData       = expData.dat( expData.dat.session > 0 , :);  % TODO: change session number
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

initCond      = getInitialConditions(formattedData);

n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now

%%
load('./Data/ModelFitQualityTest_Akash_v4.mat');
errBins = -90:3:90;
[~, idx] = min(dataToSave.result.f);
modelParams = dataToSave.result.x(idx, :);
plotFitResult_v4(formattedData, modelParams, initCond, errBins)

% opt_param_sigma_ext         = dataToSave.result.x(idx, 1:n_uncertainty_levels);
% opt_param_sigma_b           = dataToSave.result.x(idx, n_uncertainty_levels+1:2*n_uncertainty_levels);
% opt_param_scale             = dataToSave.result.x(idx, 2*n_uncertainty_levels + 1);
% opt_param_sigma_meta        = dataToSave.result.x(idx, 2*n_uncertainty_levels + 2);
% opt_param_Cc                = dataToSave.result.x(idx, 2*n_uncertainty_levels + 3);
% opt_param_guessrate         = dataToSave.result.x(idx, 2*n_uncertainty_levels + 4);
% opt_param_ori_scale         = dataToSave.result.x(idx, 2*n_uncertainty_levels + 5);
% opt_param_bias              = dataToSave.result.x(idx, 2*n_uncertainty_levels + 6);
% 
% % Display parameters
% for i =1:n_uncertainty_levels
%     fprintf("sigma_ext: Fit: %.4f C: %.3f, S: %.3f, D: %.3f\n", ...
%         opt_param_sigma_ext(i), ...
%         formattedData.uncertaintyVals(i, 1), ...
%         formattedData.uncertaintyVals(i, 2), ...
%         formattedData.uncertaintyVals(i, 3))
% end
% 
% % shape*scale = sigma_int^2
% % *opt_param_scale /sqrt(formattedData.stimulusEnergy(i))
% for i =1:n_uncertainty_levels
%     fprintf("b (sigma_int): Fit: %.4f C: %.3f, S: %.3f, D: %.3f\n", ...
%         opt_param_sigma_b(i), ...
%         formattedData.uncertaintyVals(i, 1), ...
%         formattedData.uncertaintyVals(i, 2), ...
%         formattedData.uncertaintyVals(i, 3))
% end
% 
% fprintf("Scale Fit: %.4f \n", opt_param_scale)
% fprintf("Meta Fit: %.4f \n", opt_param_sigma_meta)
% fprintf("Cc Fit: %.4f \n", opt_param_Cc)
% fprintf("GR Fit: %.4f \n", opt_param_guessrate)
% fprintf("Ori scale Fit: %.4f \n", opt_param_ori_scale)
% fprintf("Bias Amp Fit: %.4f \n", opt_param_bias)



