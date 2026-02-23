close all
clear all

addpath('C:\Users\avinash1598\Desktop\Uncertainty\ProcessModel\HumanExpDataAnalysis\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\Uncertainty\ProcessModel\LLScriptsUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\ProcessModel\PlotUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\ProcessModel\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\ProcessModel\OptimizationUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')

% expData            = load('Data\CORNFB01.mat'); % Yichao
%expData            = load('./Data/COR33.mat'); % Akash
expData            = load('./Data/COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming

fltData       = expData.dat( expData.dat.session > 0 , :); 
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

%% Cross-validate
errBins   = -90:1:90; % this is dx which might affect fitting. This value should be optimal. not too fine. not too coarse.

optParams.nStarts = 20;
optParams.hyperParamC1 = 10;
optParams.randomGuessModel = true;

cv_result = NLLCrossValidate(formattedData, errBins, 5, 6, optParams);
save('./CV_Data/cross_validation_tien_full.mat', 'cv_result');

%% NLL on test data

load('./CV_Data/cross_validation_tien_constr.mat', 'cv_result');
errBins   = -90:1:90;
% optParams.nStarts = 20;
% optParams.hyperParamC1 = 0;
% optParams.randomGuessModel = true;
nllData = computeNLL_CV(formattedData, errBins, cv_result, optParams);

% remove SD from the data and then see what happens

%
figure
subplot(2, 2, 1)
hold on
histogram(nllData.nllCovModel/100, DisplayName='cov model')
histogram(nllData.nllIndModel/100, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("Fit on test data (GT: cov model data)")

subplot(2, 2, 2)
histogram(nllData.deltaNLL) 
ylabel("Count")
xlabel("delta NLL (cov - ind)")
legend
title("Fit on test data (GT: cov model data)")

subplot(2, 2, 3)
hold on
bins = 0:30:300;
histogram(nllData.fvalsCov/100, DisplayName='cov model')
histogram(nllData.fvalsInd/100, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("fvals (train dataset)")

subplot(2, 2, 4)
hold on
bins = 0:30:300;
histogram(nllData.minfvalsCov/100, DisplayName='cov model')
histogram(nllData.minfvalsInd/100, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("bestfit fvals (train dataset)")
