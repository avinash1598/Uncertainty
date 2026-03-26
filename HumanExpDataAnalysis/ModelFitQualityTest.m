close all
clear all
% ITR 191
warning("Make sure generative model is correct. While generating simulated " + ...
    "data guess rate should only be used for low confidence reports.")

% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis/Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')

%expData            = load('./Data/COR33.mat');  % Akash
expData            = load('./Data/COR31.mat');   % Tien
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
%% Cov model

% ----------------------------------------
optParams.nStarts = 300; %300;
resultReducedCov = Optimize(formattedData, initCond, "cov", [], optParams, "reduced");

% Simulated data
[~, idx] = min(resultReducedCov.f);
params   = resultReducedCov.x(idx, :);
% params(end) = 0; % setting GR zero for now
data     = GenerateCovData( params, n_uncertainty_levels, 'reduced', ...
    formattedData.orientations', size(formattedData.theta_resp_all, 3) );
initCond = getInitialConditions(data);

optParams.nStarts = 300; %300
resultReducedCovSim = Optimize(data, initCond, "cov", [], optParams, 'reduced');

% ----------------------------------------
optParams.nStarts = 300; %300;
resultFullCov = Optimize(formattedData, initCond, "cov", [], optParams, "full");

% Simulated data
[~, idx] = min(resultFullCov.f);
params   = resultFullCov.x(idx, :);
% params(end - 2) = 0;
data     = GenerateCovData(params, n_uncertainty_levels, 'full', ...
    formattedData.orientations', size(formattedData.theta_resp_all, 3));
initCond = getInitialConditions(data);

optParams.nStarts = 300; %300
resultFullCovSim = Optimize(data, initCond, "cov", [], optParams, 'full');

% ----------------------------------------
optParams.nStarts = 300; %300;
resultFullJumboCov = Optimize(formattedData, initCond, "cov", [], optParams, "full_jumbo");

% Simulated data
[~, idx] = min(resultFullJumboCov.f);
params   = resultFullJumboCov.x(idx, :);
% params(end - n_uncertainty_levels - 1) = 0;
data     = GenerateCovData(params, n_uncertainty_levels, 'full_jumbo', ...
    formattedData.orientations', size(formattedData.theta_resp_all, 3));
initCond = getInitialConditions(data);

optParams.nStarts = 300; %300
resultFullJumboCovSim = Optimize(data, initCond, "cov", [], optParams, 'full_jumbo');

%% Ind model
% ----------------------------------------
optParams.nStarts = 300; %300;
resultReducedInd = Optimize(formattedData, initCond, "ind", [], optParams, "reduced");

% Simulated data
[~, idx] = min(resultReducedInd.f);
params   = resultReducedInd.x(idx, :);
% params(end) = 0;
data     = GenerateIndData(params, n_uncertainty_levels, 'reduced', ...
    formattedData.orientations', size(formattedData.theta_resp_all, 3));
initCond = getInitialConditions(data);

optParams.nStarts = 300; %300
resultReducedIndSim = Optimize(data, initCond, "ind", [], optParams, 'reduced');

% ----------------------------------------
optParams.nStarts = 300; %300;
resultFullInd = Optimize(formattedData, initCond, "ind", [], optParams, "full");

% Simulated data
[~, idx] = min(resultFullInd.f);
params   = resultFullInd.x(idx, :);
% params(end - 2) = 0;
data     = GenerateIndData(params, n_uncertainty_levels, 'full', ...
    formattedData.orientations', size(formattedData.theta_resp_all, 3));
initCond = getInitialConditions(data);

optParams.nStarts = 300; %300
resultFullIndSim = Optimize(data, initCond, "ind", [], optParams, 'full');

% ----------------------------------------
optParams.nStarts = 300; %300;
resultFullJumboInd = Optimize(formattedData, initCond, "ind", [], optParams, "full_jumbo");

% Simulated data
[~, idx] = min(resultFullJumboInd.f);
params   = resultFullJumboInd.x(idx, :);
% params(end - n_uncertainty_levels - 1) = 0;
data     = GenerateIndData(params, n_uncertainty_levels, 'full_jumbo', ...
    formattedData.orientations', size(formattedData.theta_resp_all, 3));
initCond = getInitialConditions(data);

optParams.nStarts = 300; %300
resultFullJumboIndSim = Optimize(data, initCond, "ind", [], optParams, 'full_jumbo');

%% Data structures to save
dataToSave.resultReducedCov   = resultReducedCov;
dataToSave.resultFullCov      = resultFullCov;
dataToSave.resultFullJumboCov = resultFullJumboCov;
dataToSave.resultReducedInd   = resultReducedInd;
dataToSave.resultFullInd      = resultFullInd;
dataToSave.resultFullJumboInd = resultFullJumboInd;

dataToSave.resultReducedCovSim   = resultReducedCovSim;
dataToSave.resultFullCovSim      = resultFullCovSim;
dataToSave.resultFullJumboCovSim = resultFullJumboCovSim;
dataToSave.resultReducedIndSim   = resultReducedIndSim;
dataToSave.resultFullIndSim      = resultFullIndSim;
dataToSave.resultFullJumboIndSim = resultFullJumboIndSim;

save('ModelFitQualityTest_Tien.mat', "dataToSave");
