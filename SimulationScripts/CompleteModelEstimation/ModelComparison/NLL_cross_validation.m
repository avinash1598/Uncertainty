close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/OptimizationUtils/')

% modelData                            = load('../modelContOriData_cov.mat');
modelData                            = load('../modelContOriData.mat');
% modelData.data.resp_err_all          = modelData.data.resp_err_all;
% modelData.data.confidence_report_all = modelData.data.confidence_report_all;

%% Cross-validate
errBins   = -90:0.1:90;
cv_result = NLLCrossValidate(modelData.data, errBins);
save('./CV_Data/cross_validation_ind_model.mat', 'cv_result'); % TODO: chnage it to ind model

%% NLL on test data
errBins   = -90:0.1:90;

% cv_result = load('./CV_Data/cross_validation_ind_model.mat');

optParams.nStarts = 20;
optParams.hyperParamC1 = 0;
optParams.randomGuessModel = true;
nllData = computeNLL_CV(modelData.data, errBins, cv_result.cv_result); 

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
title("Fit on test data")

subplot(2, 2, 2)
histogram(nllData.deltaNLL) 
ylabel("Count")
xlabel("delta NLL (cov - ind)")
legend
title("Fit on test data")

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


