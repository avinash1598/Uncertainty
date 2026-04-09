close all
clear all
restoredefaultpath

% ITR 191
warning("Make sure generative model is correct. While generating simulated " + ...
    "data guess rate should only be used for low confidence reports.")

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsNoBin\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\OptimizationScriptsNoBin\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% cov and singlystochastic for now

% expData            = load('../Data/COR33.mat');   % Akash
% expData            = load('../Data/COR31.mat');   % Tien
% expData            = load('./Data/COR32.mat');    % Jiaming
% expData            = load('../Data/CORNFB02.mat');  % Jonathan
% expData            = load('../Data/CORNFB01.mat');  % Yichao
expData            = load('../Data/CORNFB03.mat');  % David
% For subject 6 => David: make sure to change the orientation dependent
% error

fltSessionIdx = [0, 0, 0, 0, 2, 1]; % Last one is for David
fltData       = expData.dat( expData.dat.session > 1 , :);  % TODO: change session number
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

initCond             = getInitialConditions(formattedData, false); % boundaryEffect = false
n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now
%% Cov model

% % ----------------------------------------
% const scale
modelType = "cov"; %"singlyStochastic"; %"cov"; % 
optParams.nStarts = 30; %300;
resultFullCov = Optimize(formattedData, initCond, modelType, [], optParams, "full");

%%
% modelType = "cov";
% load("Yichao_cov_fit.mat", "resultFullCov");

paramsNameCov = ["sigma_s L1",...
    "sigma_s L2", ...
    "sigma_s L3", ...
    "sigma_s L4", ...
    "sigma_s L5", ...
    "sigma_s L6", ...
    "scale", "sigma_meta", "Cc", "exp", "GR", "Ori scale", "biasAmp"];

% paramsNameCov = ["sigma_s L1",...
%     "sigma_s L2", ...
%     "sigma_s L3", ...
%     "sigma_s L4", ...
%     "sigma_s L5", ...
%     "sigma_s L6", ...
%     "scale", "shape", "sigma_meta", "Cc", "GR", "Ori scale", "biasAmp"];

[~,idx] = min(resultFullCov.f);
fitParams = resultFullCov.x(idx, :);
for i=1:numel(fitParams)
    fprintf("%s, Fit: %.7f \n", paramsNameCov(i), fitParams(i))
end

% NLL computed using binned method
trlNLLs = getNLLForEachTrial(formattedData, initCond, modelType, "full", fitParams, []);
nll_binned_method = trlNLLs;

fprintf("\n\nNLL (fit) %.4f, NLL (binned method): %.4f) \n\n", ...
    resultFullCov.f(idx), sum(nll_binned_method));

%David
% cov: 8241.4146
% ind: 8249
% SS: 8268.9515

% Tien
% COV: 8796.4678, 8784.5232, 8767.9746 (GR 0.15, 1)
% IND: 8778.6861
% SS: 8808.8180 (GR 0.15), 8786.5863 (GR 1)
% Best: 8073

% Akash
% SS: 8909.6818
% IND: 8861.5184
% COV: 8851.8812, 8830.2938 (GR 1)
% Best case: 8126.6299 

% Jonathan
% COV: 8767.3606 - 8773.2721, 8759.8950 (GR 1), 8761.3211 (GR 0.15), 8763.9591 (GR 0.1)
% SS: 8772.8673, 8799.1061 (GR: 0.15), 8786.0196 (GR: 0.2), 8818.4149 (Gr 0.1),
% IND: %8772.6786 (GR 1), 8786.6121 (GR: 0.15), 8792.9594 (GR 0.1)

% Yichao
% SS: 8601.538
% IND: 8571.2607
% COV: 8571.7586 (old cov), 8566.7755 (GR: 1)

% GAMA to lognormal - try this
%%
errBins = -90:2:90;
plotFitResultSim(formattedData, fitParams, initCond, modelType, "full", errBins)
% plotFitResultSim(formattedData, fitParams, initCond, "singlyStochastic", "full", errBins)

% % Simulated data
% [~, idx] = min(resultFullCov.f);
% params   = resultFullCov.x(idx, :);
% % params(end) = 0; % setting GR zero for now
% data     = GenerateCovData( params, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3) );
% initCond = getInitialConditions(data);
% 
% optParams.nStarts = 300; %300
% resultFullCovSim = Optimize(data, initCond, "cov", [], optParams, 'full');

%% Ind model
% ----------------------------------------
% optParams.nStarts = 100; %300;
% resultFullInd = Optimize(formattedData, initCond, "ind", [], optParams, "full");

% % Simulated data
% [~, idx] = min(resultFullInd.f);
% params   = resultFullInd.x(idx, :);
% % params(end - 2) = 0;
% data     = GenerateIndData(params, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3));
% initCond = getInitialConditions(data);
% 
% optParams.nStarts = 1; %300
% resultFullIndSim = Optimize(data, initCond, "ind", [], optParams, 'full');


%% Singly Stochastic model
% optParams.nStarts = 1; %300;
% resultFullSinglyStochastic = Optimize(formattedData, initCond, "singlyStochastic", [], optParams, "full");

% SS: 8886.5879

% modelFitResults.resultFullCov              = resultFullCov;
% modelFitResults.resultFullInd              = resultFullInd;
% modelFitResults.resultFullSinglyStochastic = resultFullSinglyStochastic;
% 
% save('ModelFitQualityTest_Akash.mat', "modelFitResults");
