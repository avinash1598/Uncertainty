close all
clear all
restoredefaultpath

addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsNoBin\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\OptimizationScriptsNoBin\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')

% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% expData            = load('../Data/COR33.mat');    % Akash
% expData            = load('../Data/COR31.mat');    % Tien
% expData            = load('./Data/COR32.mat');     % Jiaming
% expData            = load('../Data/CORNFB02.mat');  % Jonathan
expData            = load('../Data/CORNFB01.mat'); % Yichao

fltData       = expData.dat( expData.dat.session > 0 , :);  % TODO: change session number
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

initCond             = getInitialConditions(formattedData);
n_uncertainty_levels = formattedData.n_uncertainty_levels;

%% Cross-validate
optParams.nStarts = 30;
K = 5; % K-fold
nPerm = 20; % no of data permutation
cv_result = NLLCrossValidate(formattedData, initCond, K, nPerm, optParams);
%%

save('./CV_Data/cross_validation_Yichao_full.mat', 'cv_result');