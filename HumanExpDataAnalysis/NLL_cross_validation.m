close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/OptimizationUtils/')

% expData            = load('Data\CORNFB01.mat'); % Yichao
expData            = load('./Data/COR33.mat'); % Akash
% expData            = load('Data\COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming

fltData       = expData.dat( expData.dat.session > 0 , :); 
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

%% Cross-validate
errBins   = -90:3:90;
cv_result = NLLCrossValidate(formattedData, errBins);
save('./CV_Data/cross_validation_akash.mat', 'Data');

%% NLL on test data
nllData = computeNLL(formattedData, errBins, cv_result);


