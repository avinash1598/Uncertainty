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
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

expData            = load('../Data/COR33.mat');     % Akash
% expData            = load('../Data/COR31.mat');   % Tien
% expData            = load('../Data/COR32.mat');   % Jiaming
% expData            = load('../Data/CORNFB02.mat');% Jonathan
% expData            = load('../Data/CORNFB01.mat');  % Yichao
% For subject 6 => David: make sure to change the orientation dependent
% error

fltSessionIdx = [0, 0, 0, 0, 2, 1];
fltData       = expData.dat( expData.dat.session > 0 , :);  % TODO: change session number
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

initCond      = getInitialConditions(formattedData);

n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now

%%
d1 = load('../Data/ModelFitQualityTest_Akash.mat');
d2 = load('../Data/ModelFitQualityTest_Akash_v3.mat');
d3 = load('../Data/ModelFitQualityTest_Akash_v4.mat');
% d1 = load('../Data/ModelFitQualityTest_Jonathan.mat');
% d2 = load('../Data/ModelFitQualityTest_Jonathan_v3.mat');
% d3 = load('../Data/ModelFitQualityTest_Jonathan_v4.mat');
% d1 = load('../Data/ModelFitQualityTest_Tien.mat');
% d2 = load('../Data/ModelFitQualityTest_Tien_v3.mat');
% d3 = load('../Data/ModelFitQualityTest_Tien_v4.mat'); % To avoid error for now
% d1 = load('../Data/ModelFitQualityTest_Yichao.mat');
% d2 = load('../Data/ModelFitQualityTest_Yichao_v4.mat');
% d3 = load('../Data/ModelFitQualityTest_Yichao_v4.mat');

plotFitQuality_v3_v4(formattedData, d1.dataToSave, d2.dataToSave, d3.dataToSave)

